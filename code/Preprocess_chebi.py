# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 14:16:59 2023

Code developed in Python version 3.9

@author: bulk007 & ozen002

Script to go from the ChEBI compound.tsv and names.tsv file, as provided by 
https://www.ebi.ac.uk/chebi/downloadsForward.do under the flat files, containing
all the compounds in ChEBI and their synonyms, to a cleaned and preprocessed
list of possible hazards in a csv format with CHEBI identifier. This version
specifically focusses on generating specific compounds and NOT groups or classes
of compounds. Takes 15-20 minutes to run.
"""

import re
import spacy
import pandas as pd # version 1.4.4

# Load ChEBI files
chebi_compounds = pd.read_csv('../data/chebi_compounds.tsv', delimiter='\t', na_filter=False)
chebi_names = pd.read_csv('../data/chebi_names.tsv', delimiter='\t', na_filter=False)

# Get the list of ChEBI compounds and decapatalize, also get the identifiers
hazard_list_compounds = chebi_compounds['NAME'].str.lower()
id_list_compounds = chebi_compounds['CHEBI_ACCESSION']
id_list_compounds = id_list_compounds.str.replace('CHEBI:','')

# Get the list with the ChEBI compound synonyms and decapatalize, also get the identifiers
hazard_list_names = chebi_names['NAME'].str.lower()
id_list_names = chebi_names['COMPOUND_ID'].astype(str)

# Merge both the list of names and identifiers
hazard_list = pd.concat([hazard_list_compounds, hazard_list_names]).reset_index(drop=True)
id_list = pd.concat([id_list_compounds,id_list_names]).reset_index(drop=True)

# Remove double and trailing white spaces
hazard_list = hazard_list.map(lambda x: re.sub(r' +', ' ', x))
hazard_list = hazard_list.map(lambda x: re.sub(r' $', '', x))

# Some entries are only mentioned with the word 'compound', 'agent', 'group' etc. 
# behind it, we can remove them, because they are classes of compounds, not 
# compounds themselves. Some words below have a space added in front or after 
# the word (atom, steroid, substituted, crown), because these words can also be
# subwords of actual compounds.
drop_list = []

irrelevant_list = ['compound', 'agent', 'drug', 'entity', 'entities', 'group',
                   'derivative' 'conjugate', 'agonist', 'antagonist', 
                   'modulator', 'pesticide', 'acaricide', 'insecticide', 
                   ' atom', 'atoms', 'molecule', 'inhibitor', 'cluster', 'anion', 
                   ' steroid', '-steroid', 'steroids', 'steroides', 
                   'fungicide', 'safener', 'fatty acid', ' lipids', 
                   'glycolipid', 'metabolite', 'catalyst', 'adjutant', 'element', 
                   'chemical', 'medication', 'primary', 'secondary', 
                   'tertiary', 'quaternary', 'biogenic', 'substituted ', 
                   'amino acid', 'crown ', 'atomic', 'contaminant', 'nutrient', 
                   'sacchar', 'residue', 'pharmaceutical', 'depressant', 
                   'mimetic', 'poison', 'hormon', 'herbicide', ' parent', 
                   'nucleo', 'lytic', 'congestant', 'psychoti', 'refrigerant', 
                   'microbicide', 'leptic','lator', 'stimulant', 'septic', 
                   'food', 'carcinogen', 'plastic', 'stabili', 'surfactant', 
                   'nutrient', 'rodenticides', 'polymer', 'explosive', 'material', 
                   'sweetener','disruptor', 'blocker', 'blocks', 
                   'blockader', 'receptor', 'orchestra', 'adjuvant', 'reductant',
                   'oxidant', 'narcotic', 'cosmetic', 'allergen', 'solution', 
                   'acceptor', 'analogue', 'radical', 'replace', 'carrier', 
                   'pathway', 'testing', 'hello', 'example', 'unknown', 
                   'mineral', 'mixture', 'thyroid', 'extract', 'ligand', 'buffer', 
                   'tracer', 'venom', ' label', 'repel', 'glass', ' alloy', 
                   ' donor', '-donor', 'fuel', 'metal ', ' metal', 'metallic',
                   'sugar', 'wurcs', 'indicator']

# If the name of the hazard is exactly one of the names listed in the irrelevant
# list with spaces around dropped, then drop the hazard on the basis of its ID
irrelevant_list_formatted = [irrelevant.strip(' ') for irrelevant in irrelevant_list]

hazard_id_concat = pd.concat([hazard_list, id_list], axis = 1)
hazard_id_concat.columns = ['Chemical Name', 'ChEBI ID']

drop_id = hazard_id_concat.loc[
    hazard_id_concat['Chemical Name'].isin(irrelevant_list_formatted), 'ChEBI ID'].unique().tolist()

hazard_id_concat_filt = hazard_id_concat.loc[
    ~hazard_id_concat['ChEBI ID'].isin(drop_id), ]

# Split the columns again
hazard_list = hazard_id_concat_filt['Chemical Name']
id_list = hazard_id_concat_filt['ChEBI ID']

# Drop other hazards which contain the one of the words in irrelevant list
# above on the basis of their names as Leonieke does
drop_list = []

for hazard in hazard_list:
    if(any(irrelevant in hazard for irrelevant in irrelevant_list) or 
       hazard.startswith('anti') or hazard.startswith('steroid')): # Remove words that start with anti and steroid, as these are also groups, e.g. antioxidant
          drop_list.append(hazard)
          
id_list = id_list[hazard_list.map(lambda x: not x in drop_list)].reset_index(drop=True)
hazard_list = hazard_list[hazard_list.map(lambda x: not x in drop_list)].reset_index(drop=True)

# Drop another round of chemicals on the basis of their IDs. Names
# looked up to drop are determined after some results were obtained for leafy greens
# and were deemed to be irrelevant or too generic
hazard_id_concat_filt = pd.concat([hazard_list, id_list], axis = 1)

drop_list_2 = ['polystyrene', 'ester', 'polyester', 'ion', 'emulsifier', 
               'biological function', 'solvent', 'essential oil', 'fertilizer', 
               'biomarker', 'dietary supplement', 'bile acid', 'virulence factor', 
               'disinfectant', 'antigen', 'urea', 'cellulose', 'organic acid']

drop_list_2.extend([entity + 's' for entity in drop_list_2])

drop_id_2 = hazard_id_concat_filt.loc[
    hazard_id_concat_filt['Chemical Name'].isin(drop_list_2), 'ChEBI ID'].unique().tolist()

hazard_id_concat_filt = hazard_id_concat_filt.loc[
    ~hazard_id_concat_filt['ChEBI ID'].isin(drop_id_2), ]
         
# Remove another round of terms on the basis of their names now - because 
# we do not expect them to be problematic on their ID-basis but the names
# allude to generic concepts because they are abbreviations or brand names
drop_list_3 = ['authority', 'home', 'homes', 'same', 'ages']

hazard_id_concat_filt = hazard_id_concat_filt.loc[
    ~hazard_id_concat_filt['Chemical Name'].isin(drop_list_3), ]

# Split the concatenated dataframe into two separate lists again for 
# the remaining operations
hazard_list = hazard_id_concat_filt['Chemical Name'].reset_index(drop = True)
id_list = hazard_id_concat_filt['ChEBI ID'].reset_index(drop = True)

# Some entries start with 'a ' or 'an ' followed by a chemical, we will remove
# the 'a ' or 'an ' from the entries
for index,hazard in enumerate(hazard_list):
    if(re.findall(r"^an?\s.*", hazard)):
        hazard_list[index] = re.sub(r"^an?\s(.+)", r"\1", hazard)
        
# Get chemicals that contain dashes in the name (e.g. 13-acetyl-deoxynivalenol) and
# also add them without the dash (e.g. 13-acetyldeoxynivalenol), accounts for
# multiple dashes by looping over the word and checking the regex again, only 
# applies to dashes connecting two words, not words and numbers
expand_list_name = []
expand_list_id = []
for hazard, identifier in zip(hazard_list, id_list):
    loop = True
    while(loop == True):
        if(re.findall(r".*[a-z]{3,}-[a-z]{3,}.*", hazard)):
            hazard = re.sub(r"(.*[a-z]{3,})-([a-z]{3,}.*)", r"\1\2", hazard)
            expand_list_name.append(hazard)
            expand_list_id.append(identifier)
        else:
            loop = False
id_list = pd.concat([id_list, pd.Series(expand_list_id)], ignore_index=True).reset_index(drop=True)     
hazard_list = pd.concat([hazard_list, pd.Series(expand_list_name)], ignore_index=True).reset_index(drop=True)

# In order to make the element and its number combination more robust,
# we add all possible combinations (e.g. polonium-210 -> polonium210, 
# polonium 210, 210-polonium, 210 polonium, polonium)
expand_list_name = []
expand_list_id = []
for hazard, identifier in zip(hazard_list, id_list):
    
    if(re.findall(r"^[a-z]+\s*\-*\d+$", hazard)):
        expand_list_name.append(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\1-\2", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\1\2", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\2-\1", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\2\1", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\1 \2", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\2 \1", hazard))
        expand_list_id.append(identifier)
        
        if(len(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\1", hazard)) > 2):
            expand_list_name.append(re.sub(r"([a-z]+)\s*\-*(\d+)", r"\1", hazard))
            expand_list_id.append(identifier)
    
    if(re.findall(r"^\d+\s*\-*[a-z]+$", hazard)):
        expand_list_name.append(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\1-\2", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\1\2", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\2-\1", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\2\1", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\1 \2", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\2 \1", hazard))
        expand_list_id.append(identifier)
        
        if(len(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\2", hazard)) > 2):
            expand_list_name.append(re.sub(r"^(\d+)\s*\-*([a-z]+)$", r"\2", hazard))
            expand_list_id.append(identifier)
            
id_list = pd.concat([id_list, pd.Series(expand_list_id)], ignore_index=True).reset_index(drop=True)     
hazard_list = pd.concat([hazard_list, pd.Series(expand_list_name)], ignore_index=True).reset_index(drop=True) 

# Make versions of compounds also more robust and add all versions (e.g. 
# mycotoxin b1 -> mycotoxin b-1, mycotoxin b 1, b1 mycotoxin, mycotoxin etc.)
expand_list_name = []
expand_list_id = []
for hazard, identifier in zip(hazard_list, id_list):
    if(re.findall(r"[a-z]+\s[a-z]{1,2}\s*\-*\d{1,2}$", hazard)):
        expand_list_name.append(re.sub(r"([a-z]+)\s([a-z]{1,2})\s*\-*(\d{1,2})$", r"\1 \2-\3", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s([a-z]{1,2})\s*\-*(\d{1,2})$", r"\1 \2 \3", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s([a-z]{1,2})\s*\-*(\d{1,2})$", r"\1 \2\3", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s([a-z]{1,2})\s*\-*(\d{1,2})$", r"\2-\3 \1", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s([a-z]{1,2})\s*\-*(\d{1,2})$", r"\2 \3 \1", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s([a-z]{1,2})\s*\-*(\d{1,2})$", r"\2\3 \1", hazard))
        expand_list_id.append(identifier)
        expand_list_name.append(re.sub(r"([a-z]+)\s([a-z]{1,2})\s*\-*(\d{1,2})$", r"\1", hazard))
        expand_list_id.append(identifier)
            
id_list = pd.concat([id_list, pd.Series(expand_list_id)], ignore_index=True).reset_index(drop=True)     
hazard_list = pd.concat([hazard_list, pd.Series(expand_list_name)], ignore_index=True).reset_index(drop=True)   

# Drop single-word entries that are not hazards, are too general or also mean 
# something else. Many of these words were chosen by checking them against
# words in the english dictionary. 
drop_list = ['all', 'has', 'can','aim', 'alls', 'bes', 'man', 'adi', 'edi', 'protein', 
              'proteins', 'impurity', 'impurities', 'toxin', 'toxins', 'water', 'h20', 
              'ltd', 'one', 'plateau', 'impose', 'beyond', 'oxygen', 'effector', 
              'beta', 'alpha', 'inorganics', 'common', 'image', 'commons', 'acid', 
              'acids', 'rest', 'light', 'light green', 'transform', 'ion', 'ionen', 'iode',
              'iones', 'did', 'dids', 'yellow', 'null', 'biological role', 'electron', 
              'cofactor', 'salt', 'neutron', 'positron', 'nucleus', 'nucleon', 
              'steroid hormone', 'short-chain fatty aldehyde', 'psychedelics',
              'amine', 'atom', 'steroid', 'hold', 'access', 'deep', 'anthelmintic',
              'anthelmintics', 'unclassifieds', 'application', 'bactericide', 
              'bactericides', 'anaesthetics', 'anesthetics', 'commotional', 'antimetabolite',
              'propellants', 'emulsifiers', 'emulgents', 'emulgent', 'vulnerary',
              'astringent', 'diuretics', 'medicament', 'farmaco', 'neurotoxin', 
              'fertiliser', 'fertilisers', 'avicides', 'purgatives', 'purgative',
              'aperients', 'aperient', 'megaphone', 'fragrance', 'endocrine', 'adrenergics',
              'anxiolytics', 'ataractics', 'milestone', 'ballistic', 'humectants',
              'vitamin', 'vitamins', 'plexiglas', 'styrofoam', 'reducer', 'reducers',
              'oxidizer', 'oxidizers', 'oxidiser', 'oxidisers', 'racemates', 'preserval',
              'defoamer', 'defoamers', 'charcoal', 'pigment', 'pigments', 'depigmentor',
              'depigmentors', 'negatron', 'alcohol', 'alcohols', 'proclaim', 'filler',
              'fillers', 'pressor', 'pressors', 'stampede', 'velocity', 'castaway',
              'deadline', 'defender', 'clearcast', 'raptor', 'beyond', 'vulture',
              'prestige', 'prestage', 'trigger', 'mini-pill', 'minipill', 'reposal',
              'essence', 'perfume', 'parfum', 'scent', 'aroma', 'arome', 'android',
              'asphalt', 'clipper', 'roundup', 'verdict', 'stipend', 'scepter',
              'stirrup', 'prophecy', 'prophecies', 'relaxin', 'retinal', 'ionomer',
              'ionomers', 'voltage', 'boltage', 'counter', 'balance', 'divinyl',
              'bivinyl', 'vinyl', 'steward', 'pegasus', 'pectin', 'pectins', 'ionones',
              'prosper', 'impulse', 'vulvate', 'aminate', 'tenuate', 'celsius',
              'proton', 'sandal', 'torque', 'formal', 'letter', 'blazer', 'corona',
              'patrol', 'condor', 'action', 'assert', 'dagger', 'empire', 'merlin',
              'manage', 'muster', 'autumn', 'tartar', 'squill', 'spray-tox', 'serval',
              'gemini', 'cypher', 'factor', 'patrol', 'aurora', 'cohort', 'reflex',
              'parlay', 'aplace', 'versed', 'equity', 'stevia', 'finish', 'stench',
              'talbot', 'antara', 'tempo', 'probe', 'gamma','arena', 'terra',
              'theta', 'glean', 'rogue', 'dozer', 'clove', 'anana', 'green', 
              'greens', 'double green', 'basic blue', 'medic', 'boson', 'quark',
              'pipe', 'pipes', 'lipid', 'lipids', 'base', 'bases', 'basen', 'amino',
              'epoxy', 'sales', 'cuban', 'flash', 'henna', 'homes', 'rally', 'midas',
              'salsa', 'cinch', 'nylon', 'hello', 'lilly', 'lacto', 'ring assembly',
              'ring assemblies', 'ring', 'rings', 'role', 'roles', 'papa', 'test', 
              'tests', 'dump', 'camp', 'camps', 'damp', 'nonmetal', 'metal',
              'metals', 'male', 'lime', 'muse', 'peep', 'peek', 'chop', 'pope', 'mope',
              'mold', 'nape', 'dean', 'fame', 'cape', 'snap', 'nova', 'decaps', 'aura',
              'gulf', 'tram', 'type', 'leap', 'sits', 'lost', 'dyes', 'tech', 'aqua',
              'perk', 'avid', 'bore', 'pete', 'stam', 'chic', 'rump', 'mess', 'sage',
              'bold', 'keto', 'chap', 'fat', 'wax', 'cpu', 'him', 'dec', 'mon', 'ash',
              'pee', 'pea', 'pan', 'mad', 'org', 'int', 'ski']
id_list = id_list[hazard_list.map(lambda x: not x in drop_list)].reset_index(drop=True)
hazard_list = hazard_list[hazard_list.map(lambda x: not x in drop_list)].reset_index(drop=True)

# Dropped the entries of length smaller than 4
id_list = id_list[hazard_list.map(lambda x: len(x) >= 4)].reset_index(drop=True)
hazard_list = hazard_list[hazard_list.map(lambda x: len(x) >= 4)].reset_index(drop=True)

# Remove entries that are only digits
id_list = id_list[hazard_list.map(lambda x: not(x.isdigit()))].reset_index(drop=True)
hazard_list = hazard_list[hazard_list.map(lambda x: not(x.isdigit()))].reset_index(drop=True)

# Make sure to add both the plural and singular forms, and since they are compounds
# we will just use the added 's' at the end for plural
expand_list_name = []
expand_list_id = []
for hazard, identifier in zip(hazard_list, id_list):
    version_match = re.search(r"\s[a-z]{1,2}\s*\-*\d{1,2}$", hazard) # check if a compound ends in a version, because then the compound needs to made plural
    if(version_match):
      
        if(hazard[:version_match.span()[0]][-1] == 's' and len(hazard[:version_match.span()[0]]) > 4):
            expand_list_name.append(hazard[:version_match.span()[0]][:-1] + hazard[version_match.span()[0]:])
            expand_list_id.append(identifier) 
        else:
            expand_list_name.append(hazard[:version_match.span()[0]] + 's' + hazard[version_match.span()[0]:])
            expand_list_id.append(identifier) 
            
    else:
        
        if(hazard[-1] == 's' and len(hazard) > 4 and not(hazard[-2:] == ' s') and not re.findall(r"\d", hazard)):
              expand_list_name.append(hazard[:-1])
              expand_list_id.append(identifier) 
        elif(not(hazard[-1] == 's') and len(hazard) > 4):
            expand_list_name.append(hazard + 's')
            expand_list_id.append(identifier) 
            
id_list = pd.concat([id_list, pd.Series(expand_list_id)], ignore_index=True).reset_index(drop=True)     
hazard_list = pd.concat([hazard_list, pd.Series(expand_list_name)], ignore_index=True).reset_index(drop=True)   

# Remove empty string
id_list = id_list[hazard_list.map(lambda x: not(x == ''))].reset_index(drop=True)
hazard_list = hazard_list[hazard_list.map(lambda x: not(x == ''))].reset_index(drop=True)

# Transform back to dataframe
processed_chebi = pd.DataFrame(zip(hazard_list, id_list), columns=['NAME', 'COMPOUND_ID']).astype('str')

# Drop the chemical names that end with 'rna' on the basis of their
# IDs - because they are RNAs
rna_pattern = re.compile(r'^.*rna$')
drop_id = list(set(processed_chebi.loc[processed_chebi.NAME.map(
    lambda row: bool(re.search(rna_pattern, row))), 'COMPOUND_ID'].tolist()))

processed_chebi = processed_chebi.loc[~processed_chebi['COMPOUND_ID'].isin(drop_id)]

# Drop rows with duplicate combination of name and id
processed_chebi = processed_chebi.drop_duplicates()

# Remove double entries of names with different ids, keep the first occurence
processed_chebi = processed_chebi.drop_duplicates(subset=['NAME'], keep='first')

# Sort the list from longest string to smallest
processed_chebi = processed_chebi.sort_values(by='NAME', key=lambda x: x.str.len(), ascending=False)

# Write to CSV
save_path = '../data/hazards_preprocessed.csv'
processed_chebi.to_csv(save_path, header=False, index=False)
