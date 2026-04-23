#!/usr/bin/env python3
#
# (c) INCOME Hackathon 2018, Bernried, Daniel^2
#

import sys
import numpy as np
import json

try: 
  import amici.sbml_import
except:
  from amici import sbml_import
  
def symengineMatrixToNumpy(x, astype='float'):
    return np.array(x).reshape(x.shape).astype(astype)

def getModelJSON(sbml_file_name):

    importer = amici.sbml_import.SbmlImporter(sbml_file_name, check_validity=False)


    observables = amici.sbml_import.assignmentRules2observables(importer.sbml,
                                                filter_function=lambda variableId:
    variableId.getId().startswith('observable_') and not
    variableId.getId().endswith('_sigma'))
    importer.processSBML()
    # importer.computeModelEquations()

    # Extract compartment information directly from libsbml.
    # We emit each compartment's ID, numeric size (if set), and spatial dimensions,
    # plus a species->compartment mapping. dMod uses these to build first-class
    # `compartments` and `compartmentOf` slots on the imported eqnlist.
    sbml = importer.sbml
    compartments = []
    for i in range(sbml.getNumCompartments()):
        c = sbml.getCompartment(i)
        compartments.append({
            'id': c.getId(),
            'size': c.getSize() if c.isSetSize() else None,
            'spatialDimensions': c.getSpatialDimensions() if c.isSetSpatialDimensions() else None,
        })
    speciesCompartments = {}
    for i in range(sbml.getNumSpecies()):
        sp = sbml.getSpecies(i)
        speciesCompartments[sp.getId()] = sp.getCompartment()

    S = symengineMatrixToNumpy(importer.stoichiometricMatrix)
    dataPy = {
        'S': importer.stoichiometricMatrix.tolist(),
        'v': [str(x) for x in importer.fluxVector],
        'p': importer.parameterIndex,
        'stateNames': symengineMatrixToNumpy(importer.symbols['species']['sym'], astype='str').tolist(),
        'parameterNames': symengineMatrixToNumpy(importer.symbols['parameter']['sym'], astype='str').tolist(),
        'x0': symengineMatrixToNumpy(importer.speciesInitial, astype='str').tolist(),
        "observables": observables,
        "compartments": compartments,
        "speciesCompartments": speciesCompartments,
    }
    data = json.dumps(dataPy)

    return data

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s SBML-FILE-NAME [OUTFILE]' % __file__)
        sys.exit(1)

    sbml_file_name = sys.argv[1]
    output = getModelJSON(sbml_file_name)

    if len(sys.argv) > 2:
        outfile = sys.argv[2]
        with open(outfile, "w") as f:
            f.write(output)
    else:
        print(output)
