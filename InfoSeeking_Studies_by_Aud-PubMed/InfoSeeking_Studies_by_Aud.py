# -*- coding: utf-8 -*-
"""
Created on Tues Nov  7 20:01:45 2017

@author: Dan Wendling

  ===========================================================
  |  PYTHON FOR BROWSING MESH 'PERSONS' TERMS IN PUBMED:    |
  |  Retrieve counts for selected branch of the MeSH tree.  |
  ===========================================================

This example counts PubMed studies that describe how various audiences 
seek and consume biomedical information. It aids in updating a wiki
page that has hyperlinks to pubmed.gov for each audience category.

Task: Let's say you are creating a new web resource for a particular 
audience, such as "caregivers." What do we know about the information
needs, information seeking, and information use behaviors of this audience? 
Run this report and you will find the number of pubmed.gov records for 
studies of this type, categorized by EACH AUDIENCE listed in the Persons 
branch of the Medical Subject Headings (MeSH) tree, 
https://www.ncbi.nlm.nih.gov/mesh/68009272.

The BioPython package must be installed from the Anaconda Environments tab.
(I know very little about this package; this is not an endorsement.)
More info: http://biopython.org/DIST/docs/tutorial/Tutorial.html

Note: Studies recently added to Medline will be left out, because 
this procedure only retrieves records that already have MeSH indexing 
assigned. You might be missing months or years of the latest studies.  
By setting the date range to 6 years back, you will probably be retrieving 
at least 5 years of records; in many cases, this will NOT include the 
current year. Contact me for ways to search that will somewhat solve this 
problem, with manual intervention from you.
"""


#%%

import pandas as pd
from Bio import Entrez
import time

# Always tell NCBI who you are. Create MyNCBI acct, use address here.
# https://www.ncbi.nlm.nih.gov/myncbi/
Entrez.email =  "A.N@example.gov"

'''
# If you need the list of database names
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
         print('Total records found = ' + row["Count"])

# Load several levels of the Persons branch of MeSH
audienceList = pd.read_csv('personsBranch.csv')

'''
# FYI's for occasional use
audienceList.head()
audienceList.columns
len(audienceList)

# As needed - when creating a new strategy, set a small sample to test with
'''

'''
testList = audienceList[2:10]
audienceList = testList
testList
'''

# Get row count; use audienceList or testList
rowCount =  audienceList.shape[0]

'''
# As needed, confirm you can loop over audienceList correctly:
for i in range(rowCount):
    print(audienceList.iloc[i]["MeSH"])
'''

# Core search strategy; everything outside term list. Escape characters as needed; 
coreStrat = '[Mesh] AND \"Information Seeking Behavior\"[Mesh] AND \"last 6 year\"[dp]'
coreStrat

# Print to screen how to use this info at pubmed.gov
searchAdvice = 'Search at pubmed.gov using "AudienceTermHere"' + coreStrat
print(searchAdvice)

'''
# If you want to test at pubmed.gov. Change first term at pubmed.gov if zero retrievals
showWholeStrat = '"' + testList.iloc[0]["MeSH"] + '"' + coreStrat
showWholeStrat
'''

# Do one pubmed search for each row; return each term and its record count
rowCount =  audienceList.shape[0]


# Create new list
termList = pd.DataFrame(columns=['MeSH', 'Count'])

linkStart = "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term="
linkEnd = "%5BMesh%5D+AND+%22Information+Seeking+Behavior%22%5BMesh%5D+AND+%22last+6+year%22%5Bdp%5D\">"

for i in range(rowCount):
    strat = audienceList.iloc[i]['MeSH'] + coreStrat
    search_handle = Entrez.esearch(db="pubmed", term=strat, rettype="count")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    currTerm = audienceList.iloc[i]["MeSH"]
    if int(search_results["Count"]) > 0: # Don't need to save terms with no result
        termList.loc[i] = linkStart + currTerm + linkEnd + currTerm + "</a>", search_results["Count"]
    print(audienceList.iloc[i]["MeSH"] + ' ' + search_results["Count"]) # not required; shows you progress
    time.sleep(.34) # Limit requests per sec or NCBI might block you!!
    # If you are using a validated MyNCBI account you can make 
    # 10 requests/sec, meaning, you could set the above to .1

# termList

# Turn off 50-character output limit
pd.set_option('display.max_colwidth', -1)

# Write to file in the working directory
termList.to_html(open('termList.html', 'w'), escape=False, index=False)

