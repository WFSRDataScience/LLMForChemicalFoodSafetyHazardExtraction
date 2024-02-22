# -*- coding: utf-8 -*-
"""
Created on Wed Feb 1 13:24:18 2023

@author: ozen002

Creates a dataframe containing all tokens from cleaned versions of abstracts 
along with some document, sentence ids and linguistic properties
"""

import time
import pandas as pd
from text_processing_funcs import abstract_to_tidy_Neris

#Read in articles / abstracts
articles = pd.read_csv('../data/abstracts_clean.csv')

idd = articles.doc_id.values

#Implement custom function to obtain df's about tokens in each article / abstract
token_df_per_abst = [abstract_to_tidy_Neris(clean_abs) 
                     for clean_abs in articles.clean_abstract.tolist()]

#Handle doc_id info so that we can use it alongside other token info
doc_ids = [[doc_id] * token_df.shape[0] for token_df, doc_id in zip(token_df_per_abst, idd)]
doc_ids_star = [doc_id for doc_id_list in doc_ids for doc_id in doc_id_list]

token_df = pd.concat(token_df_per_abst, ignore_index = True)
token_df['doc_id'] = doc_ids_star
token_df['doc_sent_id'] = token_df.apply(lambda x: str(
    x['doc_id']) + "_" + str(x['sent_id']), axis = 1)

token_df.to_csv('../data/token_df_object.csv', index = False)