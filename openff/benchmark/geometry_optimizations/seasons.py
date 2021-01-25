"""Benchmarking seasons compute spec registry.

"""

SEASONS = {
        "1:1": {
            "b3lyp-d3bj/dzvp": {
                 "method": "b3lyp-d3bj",
                 "basis": "dzvp",
                 "program": "psi4",
                 "spec_name": "b3lyp-d3bj/dzvp",
                 "spec_description": "Standard OpenFF optimization quantum chemistry specification.",
               },
            },
        "1:2": {
            "smirnoff99Frosst-1.1.0": {
                 "method": "smirnoff99Frosst-1.1.0",
                 "basis": "smirnoff",
                 "program": "openmm",
                 "spec_name": "smirnoff99Frosst-1.1.0",
                 "spec_description": "default smirnoff99Frosst-1.1.0 optimization spec",
               },
            "openff-1.0.0": {
                 "method": "openff-1.0.0",
                 "basis": "smirnoff",
                 "program": "openmm",
                 "spec_name": "openff-1.0.0",
                 "spec_description": "default openff-1.0.0 optimization spec",
               },
            "openff-1.1.1": {
                 "method": "openff-1.1.1",
                 "basis": "smirnoff",
                 "program": "openmm",
                 "spec_name": "openff-1.1.1",
                 "spec_description": "default openff-1.1.1 optimization spec",
               },
            "openff-1.2.1": {
                 "method": "openff-1.2.1",
                 "basis": "smirnoff",
                 "program": "openmm",
                 "spec_name": "openff-1.2.1",
                 "spec_description": "default openff-1.2.1 optimization spec",
               },
            "openff-1.3.0": {
                 "method": "openff-1.3.0",
                 "basis": "smirnoff",
                 "program": "openmm",
                 "spec_name": "openff-1.3.0",
                 "spec_description": "default openff-1.3.0 optimization spec",
               },
            "gaff-2.11": {
                 "method": "gaff-2.11",
                 "basis": "antechamber",
                 "program": "openmm",
                 "spec_name": "gaff-2.11",
                 "spec_description": "default gaff-2.11 optimization spec",
               }             
            }
          }

