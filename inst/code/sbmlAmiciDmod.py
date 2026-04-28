#!/usr/bin/env python3
"""SBML -> JSON adapter for dMod's `import_sbml()`.

Originally written against AMICI's symbolic SBML pipeline (Hackathon 2018);
that pipeline is now incompatible with current AMICI releases. Since the
R-side consumer (R/SBMLinterface.R) only needs a small subset of what AMICI
produced — stoichiometric matrix, kinetic-law strings, parameter values,
species initial expressions, compartment metadata — we read directly from
libsbml here. No AMICI dependency.

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

Limitations relative to the old AMICI pipeline:
- Function definitions are not pre-expanded. If the SBML uses
  <functionDefinition>, references survive verbatim in the formula string;
  dMod's downstream parser will likely fail on them. Workaround: expand
  function definitions in the SBML file beforehand (libsbml's
  `expandFunctionDefinitions()` converter does this in-place).
- AssignmentRules / RateRules are NOT inlined. For PEtab v1 problems used
  in test cases 0001-0006 this is a non-issue; advanced models that rely
  on them need the AMICI pipeline (or a libsbml-side rule-substitution
  pre-pass) once we add support.
"""

import sys
import json
import libsbml


def _formula(ast):
    """Convert a MathML AST to dMod-compatible formula text."""
    return libsbml.formulaToString(ast) if ast is not None else "0"


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
