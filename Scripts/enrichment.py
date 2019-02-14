#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:     Functional enrichment analysis
#
# Author:      Sanja
#
# Created:     20.05.2016
# Copyright:   (c) Sanja 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------


import _bioinformatics.obiGO as obiGO
import _bioinformatics.obiProb as obiProb

ontology = obiGO.Ontology()
#annotations = obiGO.Annotations("Homo sapiens", ontology=ontology)
annotations = obiGO.Annotations("yeast", ontology=ontology)

def get_enriched_terms(proteins, proteins_all):
    p_values = []
    functions = []
    res = annotations.get_enriched_terms(proteins, reference=proteins_all, prob=obiProb.Hypergeometric(), use_fdr=True, slims_only=True, aspect = ["P", "F", "C"])
    for go_id, (genes, p_value, ref) in res.items():
        print(ontology[go_id].name + " with p-value: %.4f " % p_value + ", ".join(genes))
        p_values += [p_value]
        functions += [ontology[go_id].name]
    return p_values, functions



proteins = ['ENSP00000259631', 'ACTN2', 'ACVR1', 'FNTA', 'GATA2', 'PML', 'ARF1', 'GGA3']
proteins_all = ['ENSP00000259631', 'MAP2K4', 'FLNC', 'MYPN', 'ACTN2', 'ACVR1', 'FNTA', 'GATA2', 'PML', 'ARF1', 'GGA3', 'ARF3', 'XRN1', 'LSM1', 'HLC3']

p_values, functions = get_enriched_terms(proteins, proteins_all)

