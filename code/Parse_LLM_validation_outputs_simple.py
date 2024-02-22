# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 11:11:46 2023

@author: ozen002

Parse responses from LLM to simple prompt for leafy greens and shellfish 
so that we can measure performance especially precision wise 
"""

import pandas as pd
import ast
import re

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
shellfish_df = pd.read_csv('../data/llm_outputs_shellfish.csv')

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

#We try parsing responses to simple prompt
for _, row in leafy_df.iterrows():
    
    simple = row['simple_prompt']

    #Check if there seems to be a dictionary pattern in the simple response
    if bool(re.search(r"```\n*\{.*\}\n*```", simple, flags=re.DOTALL)): 

        try: 
            simple_dict_in_str = re.search(
                r"```\n*(\{.*\})\n*```", simple, flags=re.DOTALL).group(1)
            simple_dict = ast.literal_eval(simple_dict_in_str)

        #Code above throwing error usually is caused by one of the two formats
        #addressed below
        except Exception as e: 
            
            #First format '{bladder cancer : arsenic}' - where nothing is string
            #in the dictionary
            
            #Comma is supposed to split key-value pairs and ':' supposed 
            #to split key and value in a given key-value pair
            comma_splits = [comma_split.split(':')
                for comma_split in simple_dict_in_str.split(',')]
            
            #If all splits are composed of two parts, go ahead
            if all([len(split) == 2 for split in comma_splits]): 
                
                for comma_split in comma_splits: 
                    
                    food_to_add = (comma_split[0].strip(' ').strip('{').
                                   strip('}').strip(' '))
                    
                    #Check if food CONTAINS one of the terms corresponding
                    #to leafies, if so add the chemical corresponding 
                    #as potentially hazardous
                    if ('leafy green' in food_to_add.lower() or 
                        'leafy greens' in food_to_add.lower() or 
                        'leafy vegetable' in food_to_add.lower() or 
                        'leafy vegetables' in food_to_add.lower()):                   
                    
                        chem_to_add = re.sub(
                            r"['\[\{]*([^\[\]\{\}']+)['\]\}]*", r'\1', 
                            comma_split[1].strip(' ').strip('{')
                            .strip('}').strip(' '))
                        
                        chemical_hazard.append(chem_to_add.lower())
                        abstracts_supporting_hazard.append(row['abstracts'])
                        
            #Second format example - inner dictionaries as values of the outer
            #dictionary with ''s missing for strings
            #'{fish: {se : [], hg: []}, tilapia: {heptachlor: [], mycotoxin: []}'
            #also aims to cover cases where the keys or closing
            #brackets are preceded by an \n
            elif bool(re.search('^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                                simple_dict_in_str, flags = re.DOTALL)): 
                
                while(bool(re.search(
                        '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                          simple_dict_in_str, flags = re.DOTALL))): 
                    
                    finding = re.search('^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                            simple_dict_in_str, flags = re.DOTALL).group(1)
                    
                    food_to_clean = finding.split('{')[0]
                    chem_to_clean = finding.split('{')[1]
                    
                    food_cleaned = (food_to_clean.strip(',').strip(' ').
                                    strip(':').strip("'").strip('\n '))

                    chem_cleaned = [chem.split(':')[0].strip(' ').strip("'").lower() 
                                    for chem in chem_to_clean.split(', ')]
                    
                    simple_dict_in_str = (simple_dict_in_str.rstrip(finding)) + '}'

                    #Check if food CONTAINS one of the terms corresponding
                    #to leafies, if so add the chemical corresponding 
                    #as potentially hazardous
                    if ('leafy green' in food_cleaned.lower() or
                        'leafy greens' in food_cleaned.lower() or
                        'leafy vegetable' in food_cleaned.lower() or 
                        'leafy vegetables' in food_cleaned.lower()): 
                    
                        chemical_hazard.extend(chem_cleaned)
                        abstracts_supporting_hazard.extend(
                            [row['abstracts']] * len(chem_cleaned))
    
            else:
                continue

        else:

           #Make sure that literal evaluation of a response is actually a dictionary 
           #and there are items in the dictionary
           if type(simple_dict) == dict and len(list(simple_dict.keys())): 
           
               #Check whether every value is a list and if so, whether every element
               #in each list conforms to the format of being a string
               confirm_value_list_bool = [type(value) == list 
                                           for value in simple_dict.values()]

               if all(confirm_value_list_bool): 
                   confirm_elements_of_value_str_dict = [
                       type(value_el) == str for value in simple_dict.values() 
                       for value_el in value] 

                   #When each criterium finally satisfied, add response from 
                   #simple prompt to foods, chemicals and abstracts
                   if all(confirm_elements_of_value_str_dict): 
                       
                        for key, value in simple_dict.items():                       

                            if ('leafy green' in key.lower() or
                                'leafy greens' in key.lower() or
                                'leafy vegetable' in key.lower() or
                                'leafy vegetables' in key.lower()):
                                    
                                chemical_hazard.extend(
                                    [val.lower() for val in value])
                                abstracts_supporting_hazard.extend(
                                    [row['abstracts']] * len(value))

                   else:
                       continue
                    
               #Maybe all values are a string instead
               #e.g. {'vegetables': 'sterigmatocystin', 'wheat': 'sterigmatocystin'}
               elif all([type(value) == str 
                         for value in simple_dict.values()]):
                   
                   for key, value in simple_dict.items():
                       
                       if ('leafy green' in key.lower() or
                           'leafy greens' in key.lower() or
                           'leafy vegetable' in key.lower() or
                           'leafy vegetables' in key.lower()):
                               
                           chemical_hazard.extend([value.lower()])
                           abstracts_supporting_hazard.extend(
                               [row['abstracts']])

               #Maybe there area again inner dictionaries as values of the outer
               #dictionary now without missing the ''s for strings 
               #Also apply some additional \n spaces
               #'{'fish': {'se' : [], 'hg': []}, 'tilapia': {'heptachlor': [], 'mycotoxin': []}'
               elif bool(re.search(
                       '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                       simple_dict_in_str, flags = re.DOTALL)):
                   
                   while(bool(re.search(
                            '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                              simple_dict_in_str, flags = re.DOTALL))): 
                        
                        finding = re.search(
                            '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                                simple_dict_in_str, flags = re.DOTALL).group(1)
                        
                        food_to_clean = finding.split('{')[0]
                        chem_to_clean = finding.split('{')[1]
                            
                        food_cleaned = (food_to_clean.strip(',').strip(' ').strip(':').strip("'").
                                        strip('\n '))
                            
                        chem_cleaned = [chem.split(':')[0].strip(' ').strip("'").lower() 
                                        for chem in chem_to_clean.split(', ')]
                        
                        simple_dict_in_str = (simple_dict_in_str.rstrip(finding)) + '}'

                        #Check if food CONTAINS one of the terms corresponding
                        #to leafies, if so add the chemical corresponding 
                        #as potentially hazardous
                        if ('leafy green' in food_cleaned.lower() or
                            'leafy greens' in food_cleaned.lower() or
                            'leafy vegetable' in food_cleaned.lower() or
                            'leafy vegetables' in food_cleaned.lower()):
                            
                            chemical_hazard.extend(chem_cleaned)
                            abstracts_supporting_hazard.extend(
                                [row['abstracts']] * len(chem_cleaned))
                            
               else:
                   continue
    
    else:
        continue

#Let's collect findings of chemical hazards for shellfish
chemical_hazard_s = []
abstracts_supporting_hazard_s = []

for _, row in shellfish_df.iterrows():
    
    simple = row['simple_prompt']
    
    #Check if there seems to be a dictionary pattern in the step_by_step response
    if bool(re.search(r"```\n*\{.*\}\n*```", simple, flags=re.DOTALL)): 
    
        try: 
            simple_dict_in_str = re.search(
                r"```\n*(\{.*\})\n*```", simple, flags=re.DOTALL).group(1)
            simple_dict = ast.literal_eval(simple_dict_in_str)
    
        #Code above throwing error usually is caused by one of the two formats
        #addressed below
        except Exception as e: 
            
            #First format '{bladder cancer : arsenic}' - where nothing is string
            #in the dictionary
            
            #Comma is supposed to split key-value pairs and ':' supposed 
            #to split key and value in a given key-value pair
            comma_splits = [comma_split.split(':')
                for comma_split in simple_dict_in_str.split(',')]
            
            #If all splits are composed of two parts, go ahead
            if all([len(split) == 2 for split in comma_splits]): 
                
                for comma_split in comma_splits: 
                    
                    food_to_add = (comma_split[0].strip(' ').strip('{').
                                   strip('}').strip(' '))
                    
                    #Check if food CONTAINS shellfish, if so add the chemical 
                    #corresponding as potentially hazardous
                    if 'shellfish' in food_to_add.lower():
                        
                        chem_to_add = re.sub(
                            r"['\[\{]*([^\[\]\{\}']+)['\]\}]*", r'\1', 
                            comma_split[1].strip(' ').strip('{')
                            .strip('}').strip(' '))
                        
                        chemical_hazard_s.append(chem_to_add.lower())
                        abstracts_supporting_hazard_s.append(row['abstracts'])
                        
            #Second format example - inner dictionaries as values of the outer
            #dictionary with ''s missing for strings
            #'{fish: {se : [], hg: []}, tilapia: {heptachlor: [], mycotoxin: []}'
            #also aims to cover cases where the keys or closing
            #brackets are preceded by an \n
            elif bool(re.search('^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                                simple_dict_in_str, flags = re.DOTALL)): 
                
                while(bool(re.search(
                        '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                          simple_dict_in_str, flags = re.DOTALL))): 
                    
                    finding = re.search('^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                            simple_dict_in_str, flags = re.DOTALL).group(1)
                    
                    food_to_clean = finding.split('{')[0]
                    chem_to_clean = finding.split('{')[1]
                    
                    food_cleaned = (food_to_clean.strip(',').strip(' ').
                                    strip(':').strip("'").strip('\n '))
    
                    chem_cleaned = [chem.split(':')[0].strip(' ').strip("'").lower() 
                                    for chem in chem_to_clean.split(', ')]
                    
                    simple_dict_in_str = (simple_dict_in_str.rstrip(finding)) + '}'
    
                    #Check if food CONTAINS one of the terms corresponding
                    #to shellfish, if so add the chemical corresponding 
                    #as potentially hazardous
                    if 'shellfish' in food_cleaned.lower():                    
                                                
                        chemical_hazard_s.extend(chem_cleaned)
                        abstracts_supporting_hazard_s.extend(
                            [row['abstracts']] * len(chem_cleaned))
    
            else:
                continue
    
        else:
    
           #Make sure that literal evaluation of a response is actually a dictionary 
           #and there are items in the dictionary
           if type(simple_dict) == dict and len(list(simple_dict.keys())): 
           
               #Check whether every value is a list and if so, whether every element
               #in each list conforms to the format of being a string
               confirm_value_list_bool = [type(value) == list 
                                          for value in simple_dict.values()]
    
               if all(confirm_value_list_bool): 
                   confirm_elements_of_value_str_dict = [
                       type(value_el) == str for value in simple_dict.values() 
                       for value_el in value] 
    
                   #When each criterium finally satisfied, add response from 
                   #simple prompt to foods, chemicals and abstracts
                   if all(confirm_elements_of_value_str_dict): 
                       
                        for key, value in simple_dict.items():                       
    
                            if 'shellfish' in key.lower(): 
    
                                chemical_hazard_s.extend(
                                    [val.lower() for val in value])
                                abstracts_supporting_hazard_s.extend(
                                    [row['abstracts']] * len(value))
    
                   else:
                       continue
                    
               #Maybe all values are a string instead
               #e.g. {'vegetables': 'sterigmatocystin', 'wheat': 'sterigmatocystin'}
               elif all([type(value) == str 
                         for value in simple_dict.values()]):
                   
                   for key, value in simple_dict.items():
                       
                       if type(key) == str and 'shellfish' in key.lower(): 
                           
                           chemical_hazard_s.extend([value.lower()])
                           abstracts_supporting_hazard_s.extend(
                               [row['abstracts']])
    
               #Maybe there area again inner dictionaries as values of the outer
               #dictionary now without missing the ''s for strings 
               #Also apply some additional \n spaces
               #'{'fish': {'se' : [], 'hg': []}, 'tilapia': {'heptachlor': [], 'mycotoxin': []}'
               elif bool(re.search(
                       '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                       simple_dict_in_str, flags = re.DOTALL)):
                   
                   while(bool(re.search(
                            '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                              simple_dict_in_str, flags = re.DOTALL))): 
                        
                        finding = re.search(
                            '^\{\n*([^\{\}]+\{\n*[^\{\}]+\n*\})+\n*\}$', 
                                simple_dict_in_str, flags = re.DOTALL).group(1)
                        
                        food_to_clean = finding.split('{')[0]
                        chem_to_clean = finding.split('{')[1]
                            
                        food_cleaned = (food_to_clean.strip(',').strip(' ').strip(':').strip("'").
                                        strip('\n '))
                            
                        chem_cleaned = [chem.split(':')[0].strip(' ').strip("'").lower() 
                                        for chem in chem_to_clean.split(', ')]
                        
                        simple_dict_in_str = (simple_dict_in_str.rstrip(finding)) + '}'
    
                        #Check if food CONTAINS one of the terms corresponding
                        #to shellfish, if so add the chemical corresponding 
                        #as potentially hazardous
                        if 'shellfish' in food_cleaned.lower():                    
                                           
                            chemical_hazard_s.extend(chem_cleaned)
                            abstracts_supporting_hazard_s.extend(
                                [row['abstracts']] * len(chem_cleaned))
                            
               else:
                   continue
    
    else:
        continue

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

#Also strip spaces preceding or trailing the responses as well

leafy_results_df['chem_returned_in_resp'] = leafy_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_resp: re.sub('(?<=\s)\(.*?\)$', '', chem_resp))
        
leafy_results_df['chem_returned_in_resp'] = leafy_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_resp: chem_resp.strip(' ').strip('\n').strip(' '))        
        
shellfish_results_df['chem_returned_in_resp'] = shellfish_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_resp: re.sub('(?<=\s)\(.*?\)$', '', chem_resp))

shellfish_results_df['chem_returned_in_resp'] = shellfish_results_df[
    'chem_returned_in_resp'].map(
        lambda chem_resp: chem_resp.strip(' ').strip('\n').strip(' '))        
        
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
leafy_grped_by.reset_index().to_csv('../results/validation_leafy_llm_simple_with_eval.csv', 
                                     index = False)
shellfish_grped_by.reset_index().to_csv('../results/validation_shellfish_llm_simple_with_eval.csv', 
                                         index = False)