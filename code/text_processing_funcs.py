# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 15:49:53 2023

@author: ozen002

This script contains functions that provide a customised spaCy tokenizer,
and create a dataframe of tokens seen in abstracts with some of their 
linguistic properties in the columns.
"""

import spacy
import pandas
from spacy.tokenizer import Tokenizer
from spacy.util import compile_infix_regex, compile_prefix_regex, compile_suffix_regex

#Introduce the multi-word chemicals that need to be introduced as a special
#case of single token to spaCy - they should contain space to be identified as 
#such
def create_chemical_other_tokens():
    
    chem_df = pandas.read_csv('../data/hazards_preprocessed_vNeris.csv', 
                              header = None)
    chem_df.columns = ['Chemical Name', 'ChEBI ID']
    
    chem_filt_list = chem_df.loc[chem_df['Chemical Name'].str.contains(' '), 
                                 'Chemical Name'].tolist()
    
    #Also add the food term "leafy green(s)" and "leafy vegetable(s)" to this list 
    #as we also need this to be tokenized as one token (leafy green/veg - to be used 
    #as a validation case)
    chem_filt_list.extend(['leafy green', 'leafy greens', 
                           'leafy vegetable', 'leafy vegetables'])
    
    return chem_filt_list

# *************
# Spacy utils 
# *************

def custom_tokenizer(nlp, prefix_re, suffix_re, infix_re):
    return Tokenizer(nlp.vocab, prefix_search = prefix_re.search,
                                suffix_search = suffix_re.search,
                                infix_finditer = infix_re.finditer,
                                token_match = nlp.tokenizer.token_match,
                                rules = nlp.Defaults.tokenizer_exceptions)

def spacy_instance(): 
    
    """
    We customize the tokenizer to meet our specific needs for this task. 
    We remove parantheses from sets of prefixes and suffixes so that 
    parantheses are not used as splits between tokens. We also remove hyphens 
    from infixes because some phrases could still be introduced with hyphens. 
    We also add the multi-word chemicals that we want the tokenizer to 
    recognize as single token. 
    """
    
    sp_inst = spacy.load("en_core_web_sm")
    
    #Create the chem_filt_list that we will use for specifying multi-word 
    #tokens
    chem_filt_list = create_chemical_other_tokens()
    
    #Remove the hyphen-between-letters pattern from infix patterns
    inf = list(sp_inst.Defaults.infixes)
    inf = [x for x in inf if '-|–|—|--|---|——|~' not in x] 
    infix_re = compile_infix_regex(tuple(inf))

    #Remove parantheses from prefixes and suffixes
    prefixes = sp_inst.Defaults.prefixes
    prefixes.remove(r'\(')
    prefixes.remove(r'\)')
    prefix_re = compile_prefix_regex(prefixes)

    suffixes = sp_inst.Defaults.suffixes
    suffixes.remove(r'\(')
    suffixes.remove(r'\)')
    suffix_re = compile_suffix_regex(suffixes)
    
    sp_inst.tokenizer = custom_tokenizer(sp_inst, prefix_re, suffix_re, infix_re)
    
    #Add multi-word chemicals that we want spaCy to recognize as a single token
    for chem_2_add in chem_filt_list: 
        sp_inst.tokenizer.add_special_case(chem_2_add, [{'ORTH': chem_2_add}])
    
    return (sp_inst)

nlp = spacy_instance()

# *************
# Identifying sentences and tokens
# *************

def split_sentences_Neris(doc):
    
    sent_list = [sent for sent in doc.sents]
    return sent_list

def extract_tokens_plus_meta_Neris(doc):
    '''
    Extract tokens and metadata from individual spaCy doc.
    ref. https://towardsdatascience.com/structured-natural-language-processing-with-pandas-and-spacy-7089e66d2b10
    Originally written by Gianluca Frasso - Neris removed some columns 
    which would not be necessary
    '''
    return [
        (i.i, i.text.lower(), i.lemma_, i.tag_, 
        i.pos_, i.is_stop, i.is_alpha, 
        i.is_digit, i.is_punct) for i in doc
    ]

def tidy_tokens_Neris(docs):
    
    '''
    Creates a dataframe containing all tokens and their properties from 
    a given list of sentences.
    Originally written by Gianluca Frasso - Neris removed some columns 
    which would not be necessary. 
    '''
    
    cols = [
        "sent_id", "token_order",  "token", "lemma",
        "tag", "pos", "is_stop", "is_alpha", "is_digit", "is_punct"
    ]
    
    meta_df = []
    for ix, doc in enumerate(docs):
        meta = extract_tokens_plus_meta_Neris(doc)
        meta = pandas.DataFrame(meta)
        meta.columns = cols[1:]
        meta = meta.assign(sent_id = ix).loc[:, cols]
        meta_df.append(meta)
    
    out_df = pandas.concat(meta_df) 
    # out_df["doc_id"] = doc_id
    return out_df    

def abstract_to_tidy_Neris(text): 
    
    '''
    Applies sentence splitter on a given abstract and provides the sentences
    to the function that makes a dataframe of tokens and their properties.
    Originally written by Gianluca Frasso
    '''
    
    doc = nlp(text)
    sents = [split_sentences_Neris(doc)]
    
    #Neris' note to self: this looks confusing because sents object above
    #is a list of list where inner list contains sentences from an abstract
    #But tidy_tokens_Neris is meant to go through the list of sentences it receives
    #as input separately and spit out exactly one dataframe per abstract
    
    tidy_list = [tidy_tokens_Neris(sent_i) for sent_i in sents]
    tidy_dat = pandas.concat(tidy_list)

    return(tidy_dat)
