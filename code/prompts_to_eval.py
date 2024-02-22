# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 16:29:03 2023

@author: ozen002

I have written three different prompts for identification of chemical hazards 
in foods via LLMs. The last one is inspired by the paper 
https://arxiv.org/pdf/2305.11790.pdf (Mishra et al.)
"""

simple_prompt = """
Extract foods and chemicals that are said to mentioned to be a food safety hazard for one of the foods, to contaminate one of the foods or 
to have the potential to pose risk for human health via consumption of one of the foods in the text delimited by triple backticks. 
Provide the output in dictionary format with each different food as a separate key and the chemicals that are expressed to be hazardous, 
contaminant or to have the potential to pose risk for human health via consumption of the food as the value of the key.  

I want to warn you against 
some pitfalls. First, make sure that you bring each chemical name that is mentioned 
to be a contaminant, hazard, potentially harmful for human health via food consumption, especially make sure not to skip the specific compound names. 
Another thing is if chemicals are mentioned both with their names and abbreviations, make sure to return the full name of the chemical
instead of its abbreviation. Next warning - only provide foods and chemicals that are mentioned in the text provided, do not return 
any food or chemical that is not mentioned in the text. Also, do not try to provide more specific foods or chemicals 
if the foods or chemicals in the text are only mentioned in their general category. Another thing - refrain from providing irrelevant noun 
phrases or sentences in values just because they contain chemical names, limit the values of your dictionary to the names of 
relevant chemicals. Also, return an empty dictionary if you do not identify any chemical in food as contaminant, 
hazardous, potentially harmful for human health through consumption of food. Finally, limit your answer to the dictionary, 
no other explanation or justification is necessary. 

```{}```
"""

step_by_step_prompt = """
Your task is to perform the following actions:
1. Identify the chemicals mentioned in the text below provided between triple backticks 
and collect them in a list
2. Identify the foods mentioned in the the text below provided between triple backticks 
and collect them in a list
3. Create all combinations of foods and chemicals as tuples and collect the tuples in a list
4. Go over each food-chemical combination in the list created at step 3 and look whether the 
chemical is mentioned to be a food safety hazard for that food, to contaminate that food 
or to have the potential to pose risk for human health via consumption of that food. Store each food-chemical 
pair where chemical is said to be hazardous, contaminant or to have the potential to pose risk for human health 
via consumption of the food, in a dictionary where foods are keys and chemicals that are expressed to be hazardous, 
contaminant or to have the potential to pose risk for human health via consumption of the food are values.
5. Once you go over every food-chemical pair, return the dictionary you obtained.

I want to warn you against 
some pitfalls. First, make sure that you bring each chemical name that is mentioned 
to be a contaminant, hazard, potentially harmful for human health via food consumption, especially make sure not to skip the specific compound names. 
Another thing is if chemicals are mentioned both with their names and abbreviations, make sure to return the full name of the chemical
instead of its abbreviation. Next warning - only provide foods and chemicals that are mentioned in the text provided, do not return 
any food or chemical that is not mentioned in the text. Also, do not try to provide more specific foods or chemicals 
if the foods or chemicals in the text are only mentioned in their general category. Another thing - refrain from providing irrelevant noun 
phrases or sentences in values just because they contain chemical names, limit the values of your dictionary to the names of 
relevant chemicals. Also, return an empty dictionary if you do not identify any chemical in food as contaminant, 
hazardous, potentially harmful for human health through consumption of food. Finally, limit your answer to the dictionary, 
no other explanation or justification is necessary.

Use the following format: 
Chemicals: <chemicals you identified in the text below between triple backticks>
Foods: <foods you identified in the text below between triple backticks>
Dictionary: <dictionary storing food-chemical pairs where foods are keys and 
chemicals expressed to be hazardous, contaminant or have the potential to pose risk for human health via consumption of the food item are values>

```{}```
"""

pseudocode_prompt = """
def identify_safety_hazards(text: str) -> dict:

    'Identify the chemical and food items mentioned in the provided abstract. For each combination of chemical substances and foods, look whether the chemical substance is mentioned to be a food safety hazard for that food, to contaminate that food, has possibility to pose risk for that food in future or has the potential to pose risk for human health via food chain. Keep food-chemical substance pairs where chemical substance is said to be hazardous, contaminant, potential future risk for the food or to pose risk for human health, in a dictionary where food items are keys and the hazardous, contaminant, potentially risky or potentially harmful for human health chemical substances are values. Once you go over every food-chemical pair, return the dictionary you obtained. I want to warn you against some pitfalls. First, make sure you bring each chemical substance name that is mentioned to be a contaminant, hazard, potential risk or harmful for human health, especially make sure not to skip the specific compound names. Another thing is if chemical substances are mentioned both with their names and abbreviations, make sure return the full name of the chemical instead of its abbreviation. Next warning - only provide foods and chemical subtances that are mentioned in the text provided, do not return any food or chemical substance that is not mentioned in the text. Also, do not try to provide more specific foods or chemical substances if the foods or chemicals in the text are only mentioned in their general category. Another thing - refrain from providing irrelevant noun phrases or sentences in values just because they contain chemical substance names, limit the values of your dictionary to the names of relevant chemical substances. Also, return an empty dictionary if you do not identify any chemical substance in food as contaminant, hazardous, potentially risky or harmful for human health. Finally, limit your answer to the dictionary, no other explanation or justification is necessary.'

    #Create an empty dictionary
    chemical_hazards_per_food = {{}} 
    
    #Identify the chemical items mentioned in the provided text and collect them in a list 
    chemical_list = identify_chemicals_in_text(text)
    
    #Identify the foods items mentioned in the provided text and collect them in a list
    food_list = identify_foods_in_text(text)

    #Create all combinations of food and chemical items as tuples and collect the tuples in a list 
    food_chemical_combinations = [(food, chemical) for food in food_list for chemical in chemical_list]
    
    #Go over each food-chemical combination and look whether the chemical is mentioned to be a food safety hazard for that food, to contaminate that food or have the potential to pose risk for human health via consumption of that food
    for food, chemical in food_chemical_combinations:
    
        #Store food-chemical pairs where chemical is said to be hazardous, contaminant or to have the potential to pose risk for human health via consumption of the food, in a dictionary where foods are keys and chemicals that are expressed to be hazardous, contaminant or to have the potential to pose risk for human health via consumption of the food are values
        if chemical_is_hazardous_food(food, chemical, text) and food in chemical_hazards_per_food: 
            chemical_hazards_per_food[food].append(chemical)
            
        elif chemical_is_hazardous_food(food, chemical, text) and food not in chemical_hazards_per_food:
            chemical_hazards_per_food[food] = [chemical]
            
        else: 
            continue

    return chemical_hazards_per_food

>>>identify_safety_hazards('{}')

"""

