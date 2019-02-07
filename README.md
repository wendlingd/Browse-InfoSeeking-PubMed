# Python - List counts for InfoSeeking studies in PubMed, by audience

## Health-medicine resource for building personas, journey maps, etc.

Browse 'information needs-seeking-use' studies in PubMed with Python-generated 
report

Task: Let's say you are creating a new web resource for a particular 
audience, such as "caregivers." What do you know about the information
needs, information seeking, and information use behaviors of this audience? 
The more you know, the more effective your communication will be. Much has 
been published.

Subject matter experts, product managers, public affairs staff, etc. in
medicine- and health-related disciplines could use an ongoing connection to
this type of research.

Through this basic Python script, the Persons branch of the Medical Subject 
Headings (MeSH) tree, https://www.ncbi.nlm.nih.gov/mesh/68009272, becomes a
useful tool for accessing this research by individual audience types. The
script provides a standing count of studies for each named audience, that you
can retrieve from pubmed.gov.

Run this bibliometric report periodically so you and your staff can have an 
uncomplicated foothold into this type of research.

The result is an HTML file with hyperlinks to pubmed.gov, for each Persons 
term that retrieves PubMed records added within the past 6 years. Because
the script only retrieves records assigned subject headings (MeSH), the newest,
unindexed records will not be retrieved.

The script searches all terms in the top three levels of the Persons branch
of the MeSH tree; below that, some terms were included and others were not.
Terms retrieving zero results are not included in the report.

MeSH changes; you may want to update the csv file here to match the pages
starting from https://www.ncbi.nlm.nih.gov/mesh/68009272. Most yearly updates
are done in the fall; more info: https://www.nlm.nih.gov/bsd/policy/yep_background.html.

## Sample report output

Number of studies by audience type. 

