# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 18:09:52 2023

@author: ozen002

Code to generate queries to search in EuropePMC, code originally written by
Leonieke vd Bulk and modified for purposes of this project
"""
# Import packages
import pandas as pd

# NERIS - Create the queries by combining chemical hazard terms with public health
# terms - these query terms come from the literature research terms that were
# used by Jen and her team in the NVWA report concerning the leafy greens
chem_hazard = ['food contamination', 'chemical pollutant*', 'chemical hazard*', 
               'contamina*', 'toxin*', 'toxic substance*', 'toxic compound*', 
               'pollutant*', 'agricultural chemical*', 'chemical compound*', 
               'chemical substance*', 'residu*']

# NERIS - We added some terms to 'bioaccumulation' to public_health too, because in good 
# number of abstracts about leafy greens Esther shared with Leonieke - 
# bioaccumulation came up 
public_health = ['public health', 'haccp', 'consumer protection', 'consumer*', 
                 'food safety', 'risk assessment*', 'risk analys*', 
                 'hazard analys*', 'human health*', 'health impact', 
                 'health risk*', 'bioaccumulation']  

# Note that %20 should be used instead of a space in these queries
queries = []
for pub_h in public_health:
    if(" " in pub_h):
        pub_h_space_repl = '"' + pub_h.replace(" ", "%20") + '"'
    
    else:
        pub_h_space_repl = pub_h
    
    for chem_h in chem_hazard: 
        if (" " in chem_h): 
            chem_h_space_repl = '"' + chem_h.replace(" ", "%20") + '"'
    
        else:
            chem_h_space_repl = chem_h
    
        queries.append("(ABSTRACT:{}%20OR%20TITLE:{}%20OR%20KW:{})%20AND%20"
                       "(ABSTRACT:{}%20OR%20TITLE:{}%20OR%20KW:{})"
                       .format(pub_h_space_repl, pub_h_space_repl, pub_h_space_repl,
                               chem_h_space_repl, chem_h_space_repl, chem_h_space_repl))
    
# Save the queries as a csv
queries_df = pd.DataFrame(queries)
queries_df.to_csv('../data/query_search_list_EuropePMC.csv', 
                  index = False, header = None)
