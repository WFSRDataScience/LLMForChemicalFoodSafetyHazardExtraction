# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 14:23:01 2023

@author: ozen002

We have created 3 kinds of prompts for comparison - first simpler one, second
with step-by-step instructions, third instructions provided in a function as 
a pseudocode like it's done in https://arxiv.org/pdf/2305.11790.pdf
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

#Replace 'leafy green', 'leafy greens', 'leafy vegetable' 
#and 'leafy vegetables' with their hyphenized versions
token_df['final_rep'] = np.where(
    token_df['final_rep'].isin(['leafy green', 'leafy greens', 
                                'leafy vegetable']), 
    token_df['final_rep'].str.replace(" ", "-"),
    token_df['final_rep'])

#Identify the abstracts containing leafy green / shellfish the way we exactly
#did during the word embedding process
leafy_green_list = ['leafy-green', 'leafy-greens', 'leafy-vegetable']
shellfish_list = ['shellfish']

leafy_doc_ids = list(set(token_df.loc[
    token_df['final_rep'].isin(leafy_green_list), 'doc_id'].tolist()))

shellfish_doc_ids = list(set(token_df.loc[
    token_df['final_rep'].isin(shellfish_list), 'doc_id'].tolist()))

#Let's bring the corresponding abstracts based on doc_id's we bring
leafy_abstracts = clean_abs.loc[clean_abs['doc_id'].isin(leafy_doc_ids), 
                                'clean_abstract'].tolist()

shellfish_abstracts = clean_abs.loc[clean_abs['doc_id'].isin(shellfish_doc_ids), 
                                    'clean_abstract'].tolist()

#Go over each prompt for leafy greens and shellfish, to collect results, so
#that we can show how each one of them does
leafy_prompt_responses = {'simple_prompt': [], 'step_by_step_prompt': [], 
                          'pseudo_prompt': []}

shellfish_prompt_responses = {'simple_prompt': [], 'step_by_step_prompt': [], 
                              'pseudo_prompt': []}

with torch.no_grad():

    for prompt, prompt_desc in zip([simple_prompt, step_by_step_prompt, 
                                    pseudocode_prompt], 
                                    ['simple_prompt', 'step_by_step_prompt', 
                                    'pseudo_prompt']): 

        for leafy_abstract in leafy_abstracts: 
        
            prompt_w_abs = prompt.format(leafy_abstract.strip('\n'))
            prompt_w_abs_temp = f'''### Instruction: {prompt_w_abs}
            ### Response:'''
            input_ids = tokenizer(prompt_w_abs_temp, return_tensors = 'pt').input_ids.cuda()
            output = model.generate(inputs = input_ids, max_new_tokens = 512, do_sample = False,
                                    num_beams = 1, repetition_penalty = 1.0)    
                        
            prompt_resp = tokenizer.decode(output[0])
            prompt_resp_2_print = re.search("(?<=### Response:).+", prompt_resp, 
                                            flags = re.DOTALL).group(0)
                        
            leafy_prompt_responses[prompt_desc].append(prompt_resp_2_print)

            del input_ids
            del output
            del prompt_resp
            del prompt_resp_2_print

            torch.cuda.empty_cache()
            
        for shellfish_abstract in shellfish_abstracts: 
        
            prompt_w_abs = prompt.format(shellfish_abstract.strip('\n'))
            prompt_w_abs_temp = f'''### Instruction: {prompt_w_abs}
            ### Response:'''
            input_ids = tokenizer(prompt_w_abs_temp, return_tensors = 'pt').input_ids.cuda()
            output = model.generate(inputs = input_ids, max_new_tokens = 512, do_sample = False,
                                    num_beams = 1, repetition_penalty = 1.0)    
                        
            prompt_resp = tokenizer.decode(output[0])
            prompt_resp_2_print = re.search("(?<=### Response:).+", prompt_resp, 
                                            flags = re.DOTALL).group(0)
                        
            shellfish_prompt_responses[prompt_desc].append(prompt_resp_2_print)

            del input_ids
            del output
            del prompt_resp
            del prompt_resp_2_print

            torch.cuda.empty_cache()

leafy_prompt_responses['abstracts'] = leafy_abstracts
shellfish_prompt_responses['abstracts'] = shellfish_abstracts

leafy_prompt_df = pd.DataFrame(leafy_prompt_responses)
shellfish_prompt_df = pd.DataFrame(shellfish_prompt_responses)

leafy_prompt_df.to_csv('../data/llm_outputs_leafy.csv', index = False)
shellfish_prompt_df.to_csv('../data/llm_outputs_shellfish.csv', index = False)
