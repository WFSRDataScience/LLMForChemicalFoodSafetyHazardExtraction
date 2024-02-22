# -*- coding: utf-8 -*-
"""
Created on Wed May 10 10:18:16 2023

@author: ozen002

This script filters out some of the tokens from the token dataframe depending 
on their linguistic characteristics, identifies abbreviations of chemicals, 
replaces these abbreviations with the chemical names, replaces chemical names
with their ids and finally maps the chemical id's or lemmatized versions of 
other tokens to a column called final_rep - which will be the representations
of all tokens that will be fed into the word embedding model. 
"""

import re
import numpy
import pandas
from string import punctuation

### Import data ###

# Import tokenized data
token_df = pandas.read_csv('../data/token_df_object.csv', keep_default_na = False)

# Read the data for chemicals
chem_df = pandas.read_csv('../data/hazards_preprocessed_vNeris.csv', 
                          keep_default_na = False)

chem_df.columns = ['Chemical Name', 'ChEBI_ID']

#Map chemicals to their ChEBI IDs
chem_name_id_dict = {chem_name: 'CHEBI:' + str(chebi_id) 
                    for chem_name, chebi_id 
                    in zip(chem_df['Chemical Name'], 
                           chem_df['ChEBI_ID'])}

### Filter unnecessary tokens ###

#Filter token_df
token_df_keep = token_df[
    ((token_df.token.str.startswith('(')) & 
      (token_df.token.str.endswith(')'))) |    #keep the tokens that start and end with a paranthesis
    (token_df.token.isin(chem_name_id_dict)) |  #keep the chemicals no matter what
    ((token_df.is_stop == False) &
      (token_df.is_punct == False) &
      (token_df.tag != "_SP"))  #also remove the stopwords, digits and spaces              
    ]

#Do further filtering out ops
#Remove single length tokens
token_df_keep = token_df_keep.loc[token_df_keep.token.str.len() > 1]

#We also do not want letter + punctuation and punctuation + letter combinations
#of length 2 (such as 'e.', 'm.', 'r =' etc.)

#Plus minus after a letter could denote an anion or cation, we do not delete
#letters followed by plus or minus
punct_wout_plus_minus = "".join(
    [punct for punct in [*punctuation] if punct not in ['+', '-']])

token_df_keep['undesirable_let_punct'] = token_df_keep['token'].map(
    lambda token: bool(re.search(rf'^[a-z][{punct_wout_plus_minus}]$', 
                  token)) or bool(re.search(rf'^[{punctuation}][a-z]$', token)))

token_df_keep = token_df_keep.loc[~token_df_keep['undesirable_let_punct']]

#Identify integers and numbers combined with some other weird punctuation 
#(number + non letter combinations which are not chemicals) and remove them
#from the token_df_keep dataframe
num_wout_letter_pat = r'^[^a-z]*\d+[^a-z]*$' #this would also capture integers

token_df_keep['lemma'] = numpy.where(
    ((token_df_keep.token.str.contains(num_wout_letter_pat) & 
    (token_df_keep.token.isin(list(chem_name_id_dict.keys())) == False))),     
    'number_percentage_token', 
    token_df_keep['lemma']
    )

token_df_keep = token_df_keep.loc[
    token_df_keep.lemma != 'number_percentage_token']

### Replace the potential chemical abbreviations in parantheses with the
### corresponding chemicals ###

#Indicate whether terms are surrounded by parantheses or not
token_df_keep['paranthesis_surround'] = numpy.where(
    (token_df_keep.token.str.startswith('(')) & 
    (token_df_keep.token.str.endswith(')')), True, False)

#Not each paranthesis term is eligible to be an abbreviation of a chemical
#Require presence of a letter and do not allow for space within the paranthesis
pattern = r"^\([^\s]*[a-z]+[^\s]*\)$"

token_df_keep['paran_term_fits_pattern'] = numpy.where(
    token_df_keep['paranthesis_surround'], 
    'check', 
    'not a paran term')

to_check = (token_df_keep.loc[token_df_keep['paranthesis_surround'], 'token']
            .map(lambda x: bool(re.search(pattern, x))))

for i, row in to_check.items(): 
    token_df_keep.loc[i, 'paran_term_fits_pattern'] = row

#Check whether the token before is a chemical if a token with paranthesis 
#matches a certain pattern - but to facilitate identification of previous token, 
#reset index of token_df_keep first
token_df_keep.reset_index(drop = True, inplace = True)

#Set a default value for specifying whether the token before is a chemical
token_df_keep['previous_chem'] = 'not a paran term/not fit to pattern/False'

#Specify the rows where you really have to check if the previous token is a 
#chemical (aka potential chemical rows)
to_check_2 = token_df_keep.loc[token_df_keep['paran_term_fits_pattern'] == True]

#Indices 1 before the pattern fitting terms in parantheses - aka indices for
#potential chemical terms
idx_check_pot_chem = numpy.array(to_check_2.index) - 1

#Reset_index for token_df_keep again so that you can aggregate the indices 
#corresponding to potential chemicals as a list
token_df_keep.reset_index(inplace = True) 

#Aggregate the tokens and list the indices they are located
pot_chems_idx_agg = (token_df_keep.loc[idx_check_pot_chem]
                      .groupby('token').agg({'index': list})
                      .reset_index())

#Specify the indices where potential chemical tokens are really present as 
#index_as_pot_chem and 1 further indices (where terms in parantheses reside) 
#as next_indices
pot_chems_idx_agg.columns = ['token', 'index_as_pot_chem']
pot_chems_idx_agg['next_indices'] = (pot_chems_idx_agg['index_as_pot_chem']
                                      .map(lambda x: list(numpy.array(x) + 1)))

#Strip the potential chemical terms from their opening parantheses if they also
#do not have a closing paranthesis in case they were part of a span in 
#an outer paranthesis e.g. (deoxynivalenol-3-glucoside or ((deoxynivalenol-3-glucoside
pot_chems_idx_agg['token'] = pot_chems_idx_agg['token'].map(
    lambda x: x.strip('(') if not x.endswith(')') else x)

#Maps potential chemicals to indices of terms in parantheses where terms in
#parantheses could be abbrevations for these chemicals. 
pot_chem_idx_par_term = {chem: ind_paran_term for chem, ind_paran_term 
                          in zip(pot_chems_idx_agg['token'], 
                                  pot_chems_idx_agg['next_indices'])}

#Check whether the word before term in parantheses is really a chemical
for token, idx_par in pot_chem_idx_par_term.items(): 
    
    if token in chem_name_id_dict: 
        token_df_keep.loc[idx_par, 'previous_chem'] = True
        
    #It could possibly be a chemical term wrapped in a pair of paranthesis
    elif token.strip(')').strip('(') in chem_name_id_dict:
        token_df_keep.loc[idx_par, 'previous_chem'] = True

    else: 
        continue
    
#For each document, identify the tokens in parantheses that indeed refer to
#chemicals and match them to the chemicals that they denote in a dictionary
grouped_tok_df_keep = (token_df_keep.groupby('doc_id')
                     [[col for col in token_df_keep.columns if col != 'doc_id']]
                     .agg(list))

grouped_tok_df_keep['token_moved_1_down'] = (grouped_tok_df_keep['token'].map(
    lambda row: ['[START]'] + row[:len(row) - 1]))

#Map the token in paranthesis to actual chemical names after you strip it from its
#parantheses because abbvs mostly occur without parantheses
grouped_tok_df_keep['dictionary_mapping'] = (grouped_tok_df_keep.apply(
    lambda row: {token[1:len(token) - 1] : token_1down
                 for token, token_1down, previous_chem 
                 in zip(row['token'], row['token_moved_1_down'], row['previous_chem'])
                 if previous_chem == True}, axis = 1))

grouped_tok_df_keep['modify_chem_abbvs'] = grouped_tok_df_keep.apply(
    lambda row : [token[1:(len(token)-1)] if previous_chem == True else token 
                  for token, previous_chem in 
                  zip(row['token'], row['previous_chem'])], axis = 1)

#Replace the terms that are abbreviations with their corresponding chemical
#names
grouped_tok_df_keep['modify_chem_abbvs'] = grouped_tok_df_keep.apply(
    lambda row : [row['dictionary_mapping'].setdefault(chem_abbv, token) 
                  for token, chem_abbv in 
                  zip(row['token'], row['modify_chem_abbvs'])], axis = 1)

#Replace the chemical names with their corresponding IDs
grouped_tok_df_keep['final_rep'] = grouped_tok_df_keep.apply(
    lambda row: [chem_name_id_dict.setdefault(chem_or_not, lemma) 
                 for chem_or_not, lemma in 
                 zip(row['modify_chem_abbvs'], row['lemma'])], axis = 1)

grouped_tok_df_keep = grouped_tok_df_keep[['modify_chem_abbvs', 'final_rep', 
                                           'previous_chem']]

tok_df_keep_explode = grouped_tok_df_keep.explode(
    ['modify_chem_abbvs', 'final_rep', 'previous_chem']).reset_index(names = 'doc_id')

#Drop rows with previous_chem = True because they are abbvs in parantheses - 
#the chemical occurs just before itself anyway
tok_df_keep_explode = tok_df_keep_explode.loc[
    ~(tok_df_keep_explode['previous_chem'] == True), ]

#Map each final representation to the collection of 'modify_chem_abbvs' in each
#doc_id because we are interested in the how chemicals are actually mentioned in text
rep_mapped_to_modchemabbvs = tok_df_keep_explode.groupby(
    ['doc_id', 'final_rep']).agg({'modify_chem_abbvs': list})

#Merge the two dataframes above
tok_df_keep_explode = tok_df_keep_explode.merge(
    rep_mapped_to_modchemabbvs, how = 'left', on = ['doc_id', 'final_rep'])

#Filter the dataframe to the necessary columns and change their names to more 
#intuitive names
tok_df_keep_explode = tok_df_keep_explode[[
    'doc_id', 'final_rep', 'modify_chem_abbvs_y']] 

tok_df_keep_explode.columns = ['doc_id', 'final_rep', 'original_rep_in_doc']

tok_df_keep_explode.to_csv('../data/token_df_object_filt_for_vecs.csv', 
                            index = False)