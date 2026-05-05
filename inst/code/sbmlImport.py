#!/usr/bin/env python3
"""SBML -> JSON adapter for dMod's `import_sbml()`. Reads SBML directly via
libsbml and emits a compact JSON payload with everything dMod needs:
stoichiometric matrix, kinetic-law strings, parameter values, species
initial expressions, compartment metadata.

Output JSON keys (consumed unchanged by R/SBMLinterface.R):
    S                  list of n_species rows; each row is n_reactions
                       stoichiometric coefficients (+products, -reactants).
    v                  list of n_reactions kinetic-law formula strings.
    p                  list of n_parameters numeric default values.
    parameterNames     list of n_parameters parameter ID strings.
    stateNames         list of n_species species ID strings.
    x0                 list of n_species initial-value expression strings.
                       Symbolic when an InitialAssignment exists; else the
                       numeric initialConcentration / initialAmount as text.
    observables        empty dict; PEtab observables.tsv is authoritative.
    compartments       list of {id, size, spatialDimensions} dicts.
    speciesCompartments  dict mapping species ID -> compartment ID.

Limitations:
- AssignmentRules are emitted in the JSON `assignmentRules` dict (LHS -> RHS
  formula text in L3 syntax). The R-side caller is expected to inline them
  into rates / inits / observables (see import_sbml() in R/SBMLinterface.R).
- RateRules (`<rateRule variable="X">`) are emitted in `rateRules` as
  {variable: rhs_formula}. The R-side adds a virtual reaction with
  stoichiometry +1 in row X, rate = rhs_formula, AFTER the kinetic-law
  volume-division step (RateRules already define dQ/dt directly, so no
  V-division). RateRules on non-species (parameter/compartment) are
  skipped with a warning.
- AlgebraicRules are NOT supported.
- Events (`<event>`) are emitted in `events` as a list of
  {triggerFormula, triggerTime, assignments=[{variable, formula}, ...]}.
  The R-side translates each (event, assignment) pair into a row of
  dMod's `eventlist` (var/time/value/root/method = "replace"). Triggers
  of the form `time >=/== T` or `geq(time, T)` resolve to a numeric/symbolic
  `time`; other triggers fall back to a `root` expression.
- FunctionDefinitions are inlined via libsbml's
  SBMLFunctionDefinitionConverter before the formulas are read out — this
  unwraps `Function_for_v_15(args...)` style calls in benchmark models like
  Zheng_PNAS2012 into the underlying expressions.
"""

import sys
import json
import re
import libsbml


def _formula(ast):
    """Convert a MathML AST to dMod-compatible formula text.

    Uses the SBML L3 formatter so MathML <power/> renders as `x^y` instead of
    `pow(x, y)` — R's `stats::D()` (used by dMod's symbolic Jacobian) has no
    derivative rule for `pow`. L3 also emits the time symbol as `time`, which
    matches dMod's convention.
    """
    return libsbml.formulaToL3String(ast) if ast is not None else "0"


def parse_sbml(sbml_file):
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(sbml_file)

    # Collect fatal libsbml parse errors (warnings are tolerated; many PEtab
    # benchmark models trigger spec-quibbles that don't affect semantics).
    fatals = []
    for i in range(doc.getNumErrors()):
        err = doc.getError(i)
        if err.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
            fatals.append("[%d] %s" % (err.getErrorId(), err.getMessage()))
    if fatals:
        raise RuntimeError("SBML parse failed:\n" + "\n".join(fatals))

    model = doc.getModel()
    if model is None:
        raise RuntimeError("No SBML model found in %s" % sbml_file)

    # Inline <functionDefinition> calls so downstream formulas only reference
    # builtins. Done in-place on the libsbml document.
    if model.getNumFunctionDefinitions() > 0:
        conv = libsbml.SBMLFunctionDefinitionConverter()
        conv.setDocument(doc)
        rc = conv.convert()
        if rc != libsbml.LIBSBML_OPERATION_SUCCESS:
            raise RuntimeError(
                "Failed to inline SBML <functionDefinition> elements (rc=%d)" % rc)
        model = doc.getModel()

    # --- compartments ---
    compartments = []
    for i in range(model.getNumCompartments()):
        c = model.getCompartment(i)
        compartments.append({
            'id': c.getId(),
            'size': c.getSize() if c.isSetSize() else None,
            'spatialDimensions': (int(c.getSpatialDimensions())
                                  if c.isSetSpatialDimensions() else None),
        })

    # --- species + initial-value resolution ---
    # Precedence: InitialAssignment > InitialConcentration > InitialAmount > 0.
    # We always emit a string so the R side gets a uniform character vector;
    # `import_sbml()` consumes both numeric-strings and symbolic expressions
    # via the parameter trafo.
    state_names = []
    species_compartments = {}
    x0 = []
    species_idx = {}
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        sid = sp.getId()
        state_names.append(sid)
        species_idx[sid] = i
        species_compartments[sid] = sp.getCompartment()

        ia = model.getInitialAssignment(sid)
        if ia is not None and ia.isSetMath():
            x0.append(_formula(ia.getMath()))
        elif sp.isSetInitialConcentration():
            x0.append(repr(float(sp.getInitialConcentration())))
        elif sp.isSetInitialAmount():
            x0.append(repr(float(sp.getInitialAmount())))
        else:
            x0.append("0")

    # --- parameters ---
    # SBML-explicit parameters first, then compartments (whose IDs appear as
    # parameter symbols in dMod's volume-divided kinetic laws). The dMod
    # R-side stores compartments separately (with $volume = compartment ID)
    # AND consumes parameter defaults from `p`/`parameterNames`; if a
    # compartment ID is referenced symbolically in a rate, the importer
    # needs a numeric default to fall back on for fixed-parameter wiring.
    parameter_names = []
    parameter_values = []
    for i in range(model.getNumParameters()):
        par = model.getParameter(i)
        parameter_names.append(par.getId())
        parameter_values.append(float(par.getValue()) if par.isSetValue() else 0.0)
    for c in compartments:
        if c['id'] not in parameter_names:
            parameter_names.append(c['id'])
            parameter_values.append(float(c['size']) if c['size'] is not None else 1.0)

    # --- reactions: stoichiometry + kinetic laws ---
    n_species = len(state_names)
    n_reactions = model.getNumReactions()

    # S[i][j] = stoichiometric coefficient of species i in reaction j.
    S = [[0.0] * n_reactions for _ in range(n_species)]
    flux_vector = []
    for j in range(n_reactions):
        rxn = model.getReaction(j)

        for k in range(rxn.getNumReactants()):
            sr = rxn.getReactant(k)
            sid = sr.getSpecies()
            if sid in species_idx:
                stoich = sr.getStoichiometry() if sr.isSetStoichiometry() else 1.0
                S[species_idx[sid]][j] -= float(stoich)

        for k in range(rxn.getNumProducts()):
            sr = rxn.getProduct(k)
            sid = sr.getSpecies()
            if sid in species_idx:
                stoich = sr.getStoichiometry() if sr.isSetStoichiometry() else 1.0
                S[species_idx[sid]][j] += float(stoich)

        kl = rxn.getKineticLaw()
        flux_vector.append(_formula(kl.getMath()) if kl is not None else "0")

    # --- rules (assignment + rate; algebraic rules are unsupported) ---
    # PEtab v1 SBML may use <assignmentRule> for time-varying inputs (e.g.
    # Boehm's BaF3_Epo). The R caller substitutes them into rates / inits /
    # observables; the assigned symbol then drops out of the parameter set.
    # <rateRule variable="X"> defines dX/dt = rhs and is mutually exclusive
    # with X being produced/consumed by reactions per the SBML spec — the R
    # caller adds it as a virtual reaction column on top of S/v.
    assignment_rules = {}
    rate_rules = {}
    for i in range(model.getNumRules()):
        rule = model.getRule(i)
        tc = rule.getTypeCode()
        if not rule.isSetMath():
            continue
        if tc == libsbml.SBML_ASSIGNMENT_RULE:
            assignment_rules[rule.getVariable()] = _formula(rule.getMath())
        elif tc == libsbml.SBML_RATE_RULE:
            rate_rules[rule.getVariable()] = _formula(rule.getMath())
        # AlgebraicRules silently skipped (would require a DAE solver).

    # --- events ---
    # SBML <event> has a <trigger> (a boolean expression) and a list of
    # <eventAssignment variable="V"> formulas. dMod's eventlist wants
    # (var, time, value, root, method); we emit one row per assignment.
    # `time` is extracted from triggers shaped like `time >= T`,
    # `time == T`, `geq(time, T)` (and friends). T may be numeric or a
    # symbolic parameter — both are forwarded as-is and resolved on the R
    # side. Other trigger shapes leave `triggerTime = None` and the R caller
    # falls back to a root expression.
    _time_pat_infix = re.compile(
        r'^\s*time\s*(>=|>|<=|<|==|!=)\s*(.+?)\s*$')
    _time_pat_prefix = re.compile(
        r'^\s*(geq|gt|leq|lt|eq|neq)\s*\(\s*time\s*,\s*(.+?)\s*\)\s*$')
    events = []
    for i in range(model.getNumEvents()):
        ev = model.getEvent(i)
        trig = ev.getTrigger()
        trigger_formula = (_formula(trig.getMath())
                           if trig is not None and trig.isSetMath() else None)
        trigger_time = None
        if trigger_formula:
            ft = trigger_formula
            m1 = _time_pat_infix.match(ft) or _time_pat_prefix.match(ft)
            if m1 is not None:
                rhs = m1.group(2).strip()
                # Numeric or symbolic — both valid in dMod's eventlist.
                try:
                    trigger_time = float(rhs)
                except ValueError:
                    trigger_time = rhs
        assignments = []
        for j in range(ev.getNumEventAssignments()):
            ea = ev.getEventAssignment(j)
            if not ea.isSetMath():
                continue
            assignments.append({
                'variable': ea.getVariable(),
                'formula':  _formula(ea.getMath()),
            })
        events.append({
            'id': ev.getId() if ev.isSetId() else 'event_%d' % i,
            'triggerFormula': trigger_formula,
            'triggerTime':    trigger_time,
            'assignments':    assignments,
        })

    return {
        'S': S,
        'v': flux_vector,
        'p': parameter_values,
        'parameterNames': parameter_names,
        'stateNames': state_names,
        'x0': x0,
        'observables': {},
        'compartments': compartments,
        'speciesCompartments': species_compartments,
        'assignmentRules': assignment_rules,
        'rateRules': rate_rules,
        'events': events,
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s SBML-FILE-NAME [OUTFILE]' % __file__)
        sys.exit(1)

    sbml_file_name = sys.argv[1]
    output = json.dumps(parse_sbml(sbml_file_name))

    if len(sys.argv) > 2:
        with open(sys.argv[2], 'w') as f:
            f.write(output)
    else:
        print(output)
