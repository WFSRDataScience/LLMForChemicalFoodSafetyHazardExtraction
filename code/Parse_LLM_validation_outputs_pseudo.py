# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 11:11:46 2023

@author: ozen002

Parse results from responses to pseudocode prompting of our LLM wrt leafy greens 
and shellfish so that we can measure performance
"""

import pandas as pd
import ast
import re

#Define a custom function to check if pseudo response has the desirable format - 
#returns a boolean datatype
def pseudo_resp_format_bool(pseudo_resp): 

    pseudo_resp_strip = pseudo_resp.removesuffix('</s>')

    try: 
        pseudo_resp_dict = ast.literal_eval(pseudo_resp_strip)

    #Code above throwing error, already indicates undesirable format
    except Exception as e: 
        return False
    
    else:
        #Check if pseudo_resp_dict fulfills format requirements - first is to 
        #be a dictionary and to be not empty
        if type(pseudo_resp_dict) == dict and len(list(pseudo_resp_dict.keys())): 
            
            #Check whether every value is a list and if so, whether every element 
            #in each list conforms to the format of being a string
            confirm_value_list_bool = [
                type(value) == list for value in pseudo_resp_dict.values()]

            if all(confirm_value_list_bool): 
                confirm_elements_of_value_str_dict = [
                    type(value_el) == str for value in pseudo_resp_dict.values() 
                    for value_el in value]

                if all(confirm_elements_of_value_str_dict): 
                    return True

                else:
                    return False

            else:
                return False
        else:
            return False

#Define another custom function to parse response to pseudocode prompt
def pseudo_resp_parse(pseudo_resp, leafy = True): 

    pseudo_resp_strip = pseudo_resp.removesuffix('</s>')
    pseudo_resp_strip_lower = pseudo_resp_strip.lower()
    pseudo_resp_dict = ast.literal_eval(pseudo_resp_strip_lower)

    hazardous_chemicals = []

    if leafy: 
        
        desired_key_indices = [i for i, dict_keys in 
                               enumerate(pseudo_resp_dict.keys()) 
                               if ('leafy green' in dict_keys) or 
                               ('leafy greens' in dict_keys) or
                               ('leafy vegetable' in dict_keys) or
                               ('leafy vegetables' in dict_keys)]
                               
        #If you have keys that CONTAIN the name(s) of leafy green food item
        #find the matching value lists and extend the hazardous chemicals
        #with those
        if len(desired_key_indices): 
            
            for idx in desired_key_indices: 
                
                hazardous_chemicals.extend(
                    list(pseudo_resp_dict.values())[idx])
                
    else:
        desired_key_indices = [i for i, dict_keys in 
                               enumerate(pseudo_resp_dict.keys()) 
                               if 'shellfish' in dict_keys]
        
        #If you have keys that fit the name(s) of shellfish food item
        #find the matching value lists and extend the hazardous chemicals
        #with those
        if len(desired_key_indices): 
            
            for idx in desired_key_indices: 
                
                hazardous_chemicals.extend(
                    list(pseudo_resp_dict.values())[idx])
    
    return hazardous_chemicals

#Define a custom function to implemented row-wise in a dataframe of
#responses so that if a chemical from a response is in abbreviation form
#and we need the longer name, we can access it by using a regex pattern. 
def from_abbv_to_longer(row): 

    possible_abbv = row['chem_returned_in_resp']
                
    #Check if there is a pattern where possible_abbv is
    #wrapped in a pair of parantheses
    if bool(re.search(r'.*\({}\)'.format(re.escape(possible_abbv)), 
                      row['supporting_abstract'].lower())): 
            
        #Bring the word(s) before the abbreviation wrapped in paranthesis
        abst_before_paran = re.search(r'.*(?=\({}\))'.format(
            re.escape(possible_abbv)), row['supporting_abstract'].lower()).group(0)
            
        #Obtain the last three words before to see what is the exact name
        #of chemical
        abst_before_paran_words = abst_before_paran.strip(" ").split(" ")
                
        abst_before_param_filt = abst_before_paran_words[
            len(abst_before_paran_words) - 6:]
                
        #Collect potential chemical names
        potential_chemical_names = []
        chemical_name = ''
                
        while len(abst_before_param_filt): 
                    
            chemical_name = abst_before_param_filt.pop() + ' ' + chemical_name
            potential_chemical_names.append(chemical_name.strip(' '))
        
        in_chebi = [chem for chem in potential_chemical_names
                    if chem in chem_name_id_dict]        
            
        #If a match to ChEBI exists
        if len(in_chebi):
        
            return in_chebi[-1]
        
        else:
            
            return possible_abbv

    #If the pattern does not exist anyway - treat it same as 
    #the case where you do not have any longer name that 
    #matches to ChEBI
    else:
        
        return possible_abbv
        

#Bring in inputs
leafy_df = pd.read_csv('../data/llm_outputs_leafy.csv')
llm_clean_abst_df = pd.read_csv('../data/abstracts_clean_for_llm.csv')
chem_df = pd.read_csv('../data/hazards_preprocessed.csv', 
                      header = None)

chem_df.columns = ['Chemical Name', 'ChEBI_ID']

#Map chemicals to their ChEBI IDs
chem_name_id_dict = {chem_name: 'CHEBI:' + str(chebi_id) 
                    for chem_name, chebi_id 
                    in zip(chem_df['Chemical Name'], 
                           chem_df['ChEBI_ID'])}

#Let's collect findings of chemical hazards for leafy greens
chemical_hazard = []
abstracts_supporting_hazard = []

#We pick the pseudo_prompt to evaluate because it is simpler to parse
#and has only somewhat a drop in precision performance
for _, row in leafy_df.iterrows():
    
    pseudo_response = row['pseudo_prompt']

    #If the format of response is as we wish, let's check to see
    #what hazards it has returned 
    if pseudo_resp_format_bool(pseudo_response): 
        
        hazards_returned = pseudo_resp_parse(pseudo_response)
        
        if len(hazards_returned): 
            
            chemical_hazard.extend(hazards_returned)
            abstracts_supporting_hazard.extend(
                [row['abstracts']] * len(hazards_returned))

#Let's collect findings of chemical hazards for shellfish
chemical_hazard_s = []
abstracts_supporting_hazard_s = []

shellfish_df = pd.read_csv('../data/llm_outputs_shellfish.csv')

for _, row in shellfish_df.iterrows():
    
    pseudo_response = row['pseudo_prompt']

    #If the format of response is as we wish, let's check to see
    #what hazards it has returned 
    if pseudo_resp_format_bool(pseudo_response): 
        
        hazards_returned = pseudo_resp_parse(pseudo_response, leafy = False)
        
        if len(hazards_returned): 
            
            chemical_hazard_s.extend(hazards_returned)
            abstracts_supporting_hazard_s.extend(
                [row['abstracts']] * len(hazards_returned))
        
#Identify the ones already in ChEBI and only identify
#other ones wrt paranthesis - but first make a dataframe       
leafy_results_df = pd.DataFrame(
    {'chem_returned_in_resp': chemical_hazard, 
     'supporting_abstract': abstracts_supporting_hazard})

shellfish_results_df = pd.DataFrame(
    {'chem_returned_in_resp': chemical_hazard_s, 
     'supporting_abstract': abstracts_supporting_hazard_s})        

#Some responses are returned with the name of chemical and
#its abbreviation in parantheses next to it - let's remove
#the abbreviation in parantheses 
#e.g. polychlorinated biphenyls (pcbs) --> polychlorinated biphenyls
leafy_results_df['chem_returned_in_resp'] = leafy_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_resp: re.sub('(?<=\s)\(.*?\)$', '', chem_resp))
        
shellfish_results_df['chem_returned_in_resp'] = shellfish_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_resp: re.sub('(?<=\s)\(.*?\)$', '', chem_resp))

#Now identify the ones that are in ChEBI - introduce an extra step
#to identify longer names of chemicals which are possibly in ChEBI
#if their initial abbreviation is not and change the chem_returned_in_resp
#and chems columns accordingly and introduce a ChEBI ID if the longer
#name corresponds to a ChEBI ID
leafy_results_df['chems'] = leafy_results_df['chem_returned_in_resp'].map(
    lambda chem_in_resp: chem_name_id_dict.get(chem_in_resp, ''))

leafy_results_df['chem_returned_in_resp'] = leafy_results_df.apply(
    lambda row: from_abbv_to_longer(row) if row['chems'] == '' 
    else (row['chem_returned_in_resp']), axis = 1)

leafy_results_df['chems'] = leafy_results_df['chem_returned_in_resp'].map(
    lambda chem_in_resp: chem_name_id_dict.get(chem_in_resp, ''))

shellfish_results_df['chems'] = shellfish_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_in_resp: chem_name_id_dict.get(chem_in_resp, ''))   
    
shellfish_results_df['chem_returned_in_resp'] = shellfish_results_df.apply(
    lambda row: from_abbv_to_longer(row) if row['chems'] == '' 
    else (row['chem_returned_in_resp']), axis = 1)        
        
shellfish_results_df['chems'] = shellfish_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_in_resp: chem_name_id_dict.get(chem_in_resp, ''))

#Obviously drop the rows that will not serve you
leafy_results_df = leafy_results_df.loc[
    leafy_results_df['chems'] != '']

shellfish_results_df = shellfish_results_df.loc[
    shellfish_results_df['chems'] != '']

#Match the abstracts to their doi's
leafy_dois = []
shellfish_dois = []

for _, row in leafy_results_df.iterrows(): 
    
    abst = row['supporting_abstract']

    if llm_clean_abst_df.loc[llm_clean_abst_df['clean_abstract'] == abst, 
                              'doi'].shape[0] == 1: 

        doi_of_interest = llm_clean_abst_df.loc[
            llm_clean_abst_df['clean_abstract'] == abst, 'doi'].item()

        leafy_dois.append(doi_of_interest)
        
    else: 
        leafy_dois.append('')
        
for _, row in shellfish_results_df.iterrows(): 
    
    abst = row['supporting_abstract']

    if llm_clean_abst_df.loc[llm_clean_abst_df['clean_abstract'] == abst, 
                              'doi'].shape[0] == 1: 

        doi_of_interest = llm_clean_abst_df.loc[
            llm_clean_abst_df['clean_abstract'] == abst, 'doi'].item()

        shellfish_dois.append(doi_of_interest)
        
    else: 
        shellfish_dois.append('')
        
#Provide the dois above as a column and drop the abstracts
leafy_results_df['dois'] = leafy_dois
shellfish_results_df['dois'] = shellfish_dois

leafy_results_df.drop('supporting_abstract', axis = 1, inplace = True)

shellfish_results_df.drop('supporting_abstract', axis = 1, inplace = True)

#Let's groupby chems column
leafy_grped_by = leafy_results_df.groupby('chems').agg({
    'chem_returned_in_resp': lambda chem_in_resp: '|'.join(chem_in_resp), 
    'dois': lambda doi: '|'.join(doi)})

shellfish_grped_by = shellfish_results_df.groupby('chems').agg({
    'chem_returned_in_resp': lambda chem_in_resp: '|'.join(chem_in_resp), 
    'dois': lambda doi: '|'.join(doi)})

#Add an empty row for manual evaluation
leafy_grped_by['eval'] = pd.Series([''] * leafy_grped_by.shape[0])
shellfish_grped_by['eval'] = pd.Series([''] * shellfish_grped_by.shape[0])

#Also save these as csv's
leafy_grped_by.reset_index().to_csv('../results/validation_leafy_llm_pseudo_with_eval.csv', 
                                    index = False)
shellfish_grped_by.reset_index().to_csv('../results/validation_shellfish_llm_pseudo_with_eval.csv', 
                                        index = False)