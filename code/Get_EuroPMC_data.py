# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 20:51:21 2023

@author: ozen002

Script to collect information on journal articles from EuroPMC based on given queries.

Code originally written by Leonieke vd Bulk. 
"""

# Import packages
import json
import time
import requests
import pandas as pd

# Set path to csv file contaning EuroPMC queries
queries = pd.read_csv('../data/query_search_list_EuropePMC.csv', header = None).to_numpy()

# Set boolean for the query to initialize new output file and initialize
# an empty dataframe to store all the data
initialize_file = True
data_df = pd.DataFrame(columns=['Query', 'DOI', 'Title', 'Abstract', 'PubYear'])

start_time = time.time()

for q in queries:
    
    # Get the query as a string and set the page as * to indicate the first page
    query = q[0]
    cMark = "*"
    last_page = False

    # Loop over the pages until it is the final one
    while(not last_page):
        
        # Send a http request to the REST webservice (ask for json format)
        uri = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={}&resultType=core&synonym=TRUE&cursorMark={}&pageSize=1000&format=json'.format(query, cMark)
        req = requests.get(uri)
        body = req.text
        
        # Parse the json file
        json_dict = json.loads(body)
        page_data = []
        
        # Loop over the articles in the page and get the DOI, title, abstract
        # and pubyear
        for result in json_dict['resultList']['result']:
            
            #NERIS - Do not accept abstracts later than April 2, 2023 for 
            #reproducibility of the results
            parse_pub_date = [int(str_number) for str_number in result[
                'firstPublicationDate'].split('-')]
            
            if parse_pub_date[0] > 2023:
                continue

            elif (parse_pub_date[0] == 2023) & (parse_pub_date[1] > 4): 
                continue

            elif (parse_pub_date[0] == 2023) & (parse_pub_date[1] == 4) & (parse_pub_date[2] > 2): 
                continue

            else: 
                result_row = []
                result_row.append(query)
                
                if('doi' in result):
                    result_row.append(result['doi'])
                else:
                    result_row.append('')
                
                if('title' in result):
                    result_row.append(result['title'])
                else:
                    result_row.append('')
                
                if('abstractText' in result):
                    result_row.append(result['abstractText'])
                else:
                    result_row.append('')
                
                if('pubYear' in result):
                    result_row.append(result['pubYear'])
                else:
                    result_row.append('')
                
                # Append the data to our list
                page_data.append(result_row)
        
        # Change the list to a pandas dataframe and compare for duplicates
        page_df = pd.DataFrame(page_data, columns=['Query', 'DOI', 'Title', 'Abstract', 'PubYear'])
        duplicates = page_df['DOI'].isin(data_df['DOI'])
        page_df.drop(page_df[duplicates].index, inplace=True)
        data_df = pd.concat([data_df, page_df], ignore_index = True)
        
        # If it is the first query, initialize file and write to it. Otherwise,
        # just append to the file.
        if(initialize_file):
            page_df.to_csv('../data/europmc_abstracts_food_safety_2023_bioaccum_incl.csv', 
                            index = False, header = False)
            initialize_file = False
        else:
            page_df.to_csv('../data/europmc_abstracts_food_safety_2023_bioaccum_incl.csv', 
                            index = False, header = False, mode = 'a') 
        
        # Check if we are on the last page by checking the nextCursorMark
        if(not 'nextCursorMark' in json_dict):
            last_page = True
        else:
            cMark = json_dict['nextCursorMark']

end_time = time.time()
print('It has taken {} minutes to run the first part'.format(
    (end_time - start_time) / 60))
