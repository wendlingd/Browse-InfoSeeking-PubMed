# -*- coding: utf-8 -*-
"""
Created on Tues Nov  7 20:01:45 2017

Last modified: 2019-02-07

@author: Dan Wendling, https://github.com/wendlingd/Browse-InfoSeeking-PubMed


-- EMAIL ADDRESS FOR A MyNCBI ACCOUNT REQUIRED. MORE INFO BELOW. --


  =========================================================
  |  PYTHON FOR BROWSING MESH 'PERSONS' TERMS IN PUBMED   |
  =========================================================

This report counts PubMed studies that describe how various audiences 
seek and consume biomedical information.

Requires:
    
- The BioPython package must be installed; http://biopython.org/DIST/docs/tutorial/Tutorial.html
- Email address from a MyNCBI account; https://www.ncbi.nlm.nih.gov/books/NBK3842/

Records in the pubmed.gov database often arrive without MeSH indexing; this
script requires MeSH, so some new records that can be found with searches at
pubmed.gov will not appear here. This search strategy uses 6 years as the 
cutoff, but this means you will probably retrieve at least 5 years of records.

MeSH changes; you may want to update the csv file here to match the pages
starting from https://www.ncbi.nlm.nih.gov/mesh/68009272. Most yearly updates
are done by January; more info: https://www.nlm.nih.gov/bsd/policy/yep_background.html.
"""


#%%
# ============================================
# Get started
# ============================================

import pandas as pd
from Bio import Entrez
import time
import os
import datetime
from datetime import datetime
from datetime import timedelta
today = datetime.today().strftime('%Y-%m-%d')

# If you want to reset the working directory
# os.chdir('/Users/username/Projects/pubmed')

# Always tell NCBI who you are. Create MyNCBI acct, use address here.
# https://www.ncbi.nlm.nih.gov/myncbi/
Entrez.email =  "yourRegisteredEmail@yoursite.org"

'''
# If you want the list of database names
handle = Entrez.einfo()
result = handle.read()
handle.close()
print(result)
'''

# Foundational search strategy (unpaired), total record count
totRecords = Entrez.egquery(term='"Information Seeking Behavior"[Mesh] AND "last 6 year"[dp]"')
record = Entrez.read(totRecords)
totRecords.close()
for row in record["eGQueryResult"]: # eGQueryResult returns the record count
     if row["DbName"]=="pubmed":
         print('Total records found = ' + row["Count"] + '\n')

# Load several levels of the Persons branch of MeSH
audienceList = pd.read_csv('personsBranch.csv')


'''
# As needed, confirm you can loop over audienceList correctly:
for i in range(rowCount):
    print(audienceList.iloc[i]["MeSH"])
'''

# Components of search strategy that are outside audienceList. Escape characters as needed; 
coreStrat = '[Mesh] AND \"Information Seeking Behavior\"[Mesh] AND \"last 6 year\"[dp]'
coreStrat

'''
# To reproduce searches at pubmed.gov, print this to screen
searchAdvice = 'Search at pubmed.gov using "AudienceTermHere"' + coreStrat
print(searchAdvice)
'''

# Get row count; use audienceList or testList
rowCount =  audienceList.shape[0]

# HTML parts
linkStart = '<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term='
linkEnd = '%5BMesh%5D+AND+%22Information+Seeking+Behavior%22%5BMesh%5D+AND+%22last+6+year%22%5Bdp%5D\">'


#%%
# ================================================
# Build dataframe by querying pubmed.gov database
# ================================================

# New df
databaseResult = pd.DataFrame(columns=['Indent', 'MeSH', 'Count'])

'''
To test with small group

audienceList = audienceList[0:9]
rowCount = 9
'''

# Process rows of audienceList
for i in range(rowCount):
    strat = audienceList.iloc[i]['MeSH'] + coreStrat
    search_handle = Entrez.esearch(db="pubmed", term=strat, rettype="count")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    currIndent = audienceList.iloc[i]["Indent"]
    currTerm = audienceList.iloc[i]["MeSH"]
    if int(search_results["Count"]) > 0: # Don't need to save terms with no result
        databaseResult.loc[i] = currIndent, currTerm, linkStart + currTerm + linkEnd + search_results["Count"] + '</a>'
    print(audienceList.iloc[i]["MeSH"] + ' ' + search_results["Count"]) # not required; shows you progress
    time.sleep(.34) # Limit requests per sec or NCBI might block you!!
    # If you are using a validated MyNCBI account you can make 
    # 10 requests/sec, meaning, you could set the above to .1

# Turn off 50-character output limit
pd.set_option('display.max_colwidth', -1)


#%%
# ============================================
# OUTPUT AS HTML LIST
# ============================================

htmlResult = databaseResult
htmlResult['Indent'] = htmlResult['Indent'].astype(str)
htmlResult['Count'] = htmlResult['Count'].astype(str)

htmlResult['newContent'] = '<li class="indent' + databaseResult['Indent'] + '"> ' + databaseResult['MeSH'] + ' (' + databaseResult['Count'] + ')</li>'

htmlResult = htmlResult[['newContent']]

htmlOutput = ['<html><head><title>Audience-specific information studies at pubmed.gov</title>\n<style> \
                body {margin-left:1em;} \
                ul {list-style-type: none; margin-left:0; padding-left:0;} \
                .indent1 {padding-left:0em;} \
                .indent2 {padding-left:2em;} \
                .indent3 {padding-left:4em;} \
                .indent4 {padding-left:6em;} \
                .indent5 {padding-left:8em;}</style></head><body>']
              
htmlOutput.append('<h1>Understand your audiences: Health-medical information-seeking-behavior \
                  studies at pubmed.gov</h1> \
                  <p><strong>{}</strong> - Use this report to understand information needs, information \
                  seeking, and information use behaviors for various audiences \
                  for biomedical information. Example: The \
                  &quot;caregivers&quot; link below runs the following search strategy at \
                  <a href="https://pubmed.gov">pubmed.gov</a>: Caregivers[Mesh] \
                  AND "Information Seeking Behavior"[Mesh] AND "last 6 year"[dp]. \
                  This report was built from the \
                  <a href="https://www.ncbi.nlm.nih.gov/mesh/68009272">Persons \
                  branch of the MeSH tree.</a> Terms with zero results are not shown. \
                  Note: New database records are built incrementally, and the \
                  newest records do not have MeSH terms - the newest records are \
                  not retrievable here. Counts will go out of date; Python code \
                  for running the report yourself is at \
                  <a href="https://github.com/wendlingd/Browse-InfoSeeking-PubMed"> \
                  https://github.com/wendlingd/Browse-InfoSeeking-PubMed.</a> \
                  </p><ul>'.format(today))

for index, row in htmlResult.iterrows():
    htmlOutput.append(row['newContent'])

# Convert to string
htmlString = '\n'.join(htmlOutput)

htmlEnd = '\n</ul>\n\n</body>\n</html>'
htmlString = htmlString + htmlEnd

# Write out
writeToF = open('InfoSeekingStudies.html', "w")
writeToF.write(str(htmlString))
writeToF.close()


'''
# Or, a flat table of database results (no indents)
databaseResult.to_html(open('databaseResult.html', 'w'), escape=False, index=False)
'''
