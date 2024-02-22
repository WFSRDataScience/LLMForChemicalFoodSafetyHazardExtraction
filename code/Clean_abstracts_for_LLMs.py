"""
Written by @ozen002. 

Do some text cleaning operations on the collected abstracts before we 
identify chemical hazards for food items via our LLM of choice
Different from cleaning process done for LLM, this is closer to
natural language in its written form.
"""

#Libraries
import re
import numpy
import pandas

def clean_text(text): 
    
    """
    We have a different text cleaning function for the LLM because it can
    handle the intricacies of natural language better than a word embedding 
    model does. (These intricacies can sometimes be even helpful)
    """
    
    #Remove HTML tags around the text 
    #Notice that we replace all tags with empty string without any space
    no_html = (text.replace("<i>", "").replace("<b>", "").
               replace("</i>", "").replace("</b>", "").replace("<sub>", "").
               replace("</sub>", "").replace("<sup>", "").replace("</sup>", ""))
    
    #Make sure that this is not a greedy capture
    no_html_2 = re.sub(r'<h4>(.*?)</h4>', r' ', no_html)
    
    #There should be single spaces between tokens. If there is a period, it
    #should not have a space on its left side. 
    single_space = re.sub(r'\.(?!\d+|\s)', r'. ', no_html_2)
    single_space2 = re.sub(r'\s\s+(?!\.)', r' ', single_space)
    single_space3 = re.sub(r'\s+\.', r'.', single_space2)
    single_space4 = re.sub(r'\.\s\s+', r'. ', single_space3)

    #Remove the copyright declarations at the end of abstract.
    wout_copyright = re.sub("this article is protected by copyright. all rights reserved.", 
                            r'', single_space4)

    #We introduce character limit that can follow the keyword copyright or
    #copyright sign because otherwise we can delete copyright and everything that 
    #follows because it is mentioned in an abstract
    wout_copyright2 = re.sub(r'copyright(.){0, 50}', r'', wout_copyright)
    wout_copyright3 = re.sub(r'(copyright)*(:|\s)*©(.){0,50}', r'', wout_copyright2).strip()
    
    return wout_copyright3

# Import data
articles = pandas.read_csv("../data/europmc_abstracts_food_safety_2023_bioaccum_incl.csv", 
                           on_bad_lines = 'warn', 
                           header = None)

# Specify column names
articles.columns = ['query', 'doi', 'title', 'abstract', 'year']
articles['doc_id'] = numpy.arange(0, articles.shape[0])

# Drop NaN entry and duplicated abstracts 
articles.dropna(subset = 'abstract', inplace = True)
articles.drop_duplicates(['abstract'], inplace = True)

# Also drop duplicates in terms of doi, make sure it does not consider NaN doi's
# as duplicates
articles = articles[(~articles['doi'].duplicated()) | (articles['doi'].isna())]

# Sort articles according to their length
articles['abstract_len'] = articles['abstract'].map(lambda x: len(str(x)))
articles_sorted = articles.sort_values(by = 'abstract_len', ascending = True)

# Remove the ones that are shorter than 60 characters - definitely
articles_sorted = articles_sorted.loc[articles_sorted.abstract_len > 60]

# Remove the ones which contain some sort of correction, amendment etc. 
articles_sorted = articles_sorted.loc[~((articles_sorted.abstract.str.lower().str.contains(
    'this corrects the article', case = False)) |
                                        (articles_sorted.abstract.str.lower().str.contains(
                                            'this retracts the article', case = False)) |
                                        (articles_sorted.abstract.str.lower().str.contains(
                                            'an amendment to this paper', case = False)))]
  
#Drop abstract_len column helping us doing the filters above
articles_sorted.drop(['abstract_len'], axis = 1, inplace = True)

# Reset index 
articles_sorted.reset_index(drop = True, inplace = True)

articles_sorted['clean_abstract'] = articles_sorted['abstract'].map(
    lambda x: clean_text(x))

articles_sorted.to_csv('../data/abstracts_clean_for_llm.csv', index = None)