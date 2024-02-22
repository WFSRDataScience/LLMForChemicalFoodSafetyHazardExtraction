# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 14:23:01 2023

@author: ozen002

We extract hazards for salmon, maize and dairy as test cases for Esther v Asselt
to evaluate
"""

from prompts_to_eval import simple_prompt, step_by_step_prompt, pseudocode_prompt
from auto_gptq import AutoGPTQForCausalLM, BaseQuantizeConfig
from transformers import AutoTokenizer, logging
import pandas as pd
import numpy as np
import argparse
import torch
import re
import os

#Set some env variable to prevent OOM errors
torch.cuda.empty_cache()
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:22'

#Set environment variable
os.environ["TOKENIZERS_PARALLELISM"] = "false"

#Instantiate the tokenizer and the specific Transformer model
model_name_or_path = "TheBloke/Nous-Hermes-13B-GPTQ"
model_basename = "nous-hermes-13b-GPTQ-4bit-128g.no-act.order"

use_triton = True

tokenizer = AutoTokenizer.from_pretrained(model_name_or_path, use_fast=True)

model = AutoGPTQForCausalLM.from_quantized(model_name_or_path,
        use_safetensors=True,
        trust_remote_code=True,
        device="cuda:0",
        use_triton=use_triton,
        inject_fused_mlp=False,
        quantize_config=None)

model.eval()

#Prevent printing spurious transformers error when using pipeline with AutoGPTQ
logging.set_verbosity(logging.CRITICAL) 

#Bring all abstracts - we will filter them out
clean_abs = pd.read_csv('../data/abstracts_clean_for_llm.csv')

#Read the token_dat_object that we filtered previously
token_df = pd.read_csv("../data/token_df_object_filt_for_vecs.csv", 
                       keep_default_na = False)

#Replace 'dairy product', 'dairy products' with their hyphenized versions
token_df['final_rep'] = np.where(
    token_df['final_rep'].isin(['dairy product', 'dairy products', 
                                'dairy food product', 'dairy food products']), 
    token_df['final_rep'].str.replace(" ", "-"),
    token_df['final_rep'])

#Identify the abstracts containing salmon, maize, dairy ((food) product(s)) 
#the way we exactly did during the word embedding process
salmon_list = ['salmon']
maize_list = ['maize', 'corn']
dairy_list = ['dairy', 'dairy-product', 'dairy-products', 
              'dairy-food-product', 'dairy-food-products']

salmon_doc_ids = list(set(token_df.loc[
    token_df['final_rep'].isin(salmon_list), 'doc_id'].tolist()))

maize_doc_ids = list(set(token_df.loc[
    token_df['final_rep'].isin(maize_list), 'doc_id'].tolist()))

dairy_doc_ids = list(set(token_df.loc[
    token_df['final_rep'].isin(dairy_list), 'doc_id'].tolist()))

#Let's bring the corresponding abstracts based on doc_id's we bring
salmon_abstracts = clean_abs.loc[
    clean_abs['doc_id'].isin(salmon_doc_ids), 'clean_abstract'].tolist()

maize_abstracts = clean_abs.loc[
    clean_abs['doc_id'].isin(maize_doc_ids), 'clean_abstract'].tolist()

dairy_abstracts = clean_abs.loc[
    clean_abs['doc_id'].isin(dairy_doc_ids), 'clean_abstract'].tolist()

#Go over 2 of the prompts for salmon, maize and dairy, to collect results, so
#that we can show how each one of them does
salmon_prompt_responses = {'step_by_step_prompt': [], 'pseudo_prompt': []}

maize_prompt_responses = {'step_by_step_prompt': [], 'pseudo_prompt': []}

dairy_prompt_responses = {'step_by_step_prompt': [], 'pseudo_prompt': []}

with torch.no_grad():

    for prompt, prompt_desc in zip([step_by_step_prompt, pseudocode_prompt], 
                                   ['step_by_step_prompt', 'pseudo_prompt']): 

        for salmon_abstract in salmon_abstracts: 
        
            prompt_w_abs = prompt.format(salmon_abstract.strip('\n'))
            prompt_w_abs_temp = f'''### Instruction: {prompt_w_abs}
            ### Response:'''
            input_ids = tokenizer(prompt_w_abs_temp, return_tensors = 'pt').input_ids.cuda()
            output = model.generate(inputs = input_ids, max_new_tokens = 512, 
                                    do_sample = False, num_beams = 1, 
                                    repetition_penalty = 1.0)    
                        
            prompt_resp = tokenizer.decode(output[0])
            prompt_resp_2_print = re.search("(?<=### Response:).+", prompt_resp, 
                                            flags = re.DOTALL).group(0)
                        
            salmon_prompt_responses[prompt_desc].append(prompt_resp_2_print)

            del input_ids
            del output
            del prompt_resp
            del prompt_resp_2_print

            torch.cuda.empty_cache()
            
        for maize_abstract in maize_abstracts: 
        
            prompt_w_abs = prompt.format(maize_abstract.strip('\n'))
            prompt_w_abs_temp = f'''### Instruction: {prompt_w_abs}
            ### Response:'''
            input_ids = tokenizer(prompt_w_abs_temp, return_tensors = 'pt').input_ids.cuda()
            output = model.generate(inputs = input_ids, max_new_tokens = 512, 
                                    do_sample = False, num_beams = 1, 
                                    repetition_penalty = 1.0)    
                        
            prompt_resp = tokenizer.decode(output[0])
            prompt_resp_2_print = re.search("(?<=### Response:).+", prompt_resp, 
                                            flags = re.DOTALL).group(0)
                        
            maize_prompt_responses[prompt_desc].append(prompt_resp_2_print)

            del input_ids
            del output
            del prompt_resp
            del prompt_resp_2_print

            torch.cuda.empty_cache()
            
        for dairy_abstract in dairy_abstracts: 
        
            prompt_w_abs = prompt.format(dairy_abstract.strip('\n'))
            prompt_w_abs_temp = f'''### Instruction: {prompt_w_abs}
            ### Response:'''
            input_ids = tokenizer(prompt_w_abs_temp, return_tensors = 'pt').input_ids.cuda()
            output = model.generate(inputs = input_ids, max_new_tokens = 512, 
                                    do_sample = False, num_beams = 1, 
                                    repetition_penalty = 1.0)    
                        
            prompt_resp = tokenizer.decode(output[0])
            prompt_resp_2_print = re.search("(?<=### Response:).+", prompt_resp, 
                                            flags = re.DOTALL).group(0)
                        
            dairy_prompt_responses[prompt_desc].append(prompt_resp_2_print)

            del input_ids
            del output
            del prompt_resp
            del prompt_resp_2_print

            torch.cuda.empty_cache()

salmon_prompt_responses['abstracts'] = salmon_abstracts
maize_prompt_responses['abstracts'] = maize_abstracts
dairy_prompt_responses['abstracts'] = dairy_abstracts

salmon_prompt_df = pd.DataFrame(salmon_prompt_responses)
maize_prompt_df = pd.DataFrame(maize_prompt_responses)
dairy_prompt_df = pd.DataFrame(dairy_prompt_responses)

salmon_prompt_df.to_csv('../data/llm_outputs_salmon.csv', index = False)
maize_prompt_df.to_csv('../data/llm_outputs_maize.csv', index = False)
dairy_prompt_df.to_csv('../data/llm_outputs_dairy.csv', index = False)
