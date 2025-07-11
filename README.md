# Python - List counts for InfoSeeking studies in PubMed, by audience

### For creators of human-centered health and medicine information products

Task: Let's say you are creating a new web resource for a particular 
audience, such as "caregivers." What do you know about the information
needs, information seeking, and information use behaviors of this audience? 
The more you know, the more effective your communication will be. Much has 
been published. Having an ongoing connection to audience research can improve 
your effectiveness.

This script searches terms in the top three levels of the Persons branch
of the Medical Subject Headings (MeSH) tree, https://www.ncbi.nlm.nih.gov/mesh/68009272. 
In the levels below the third, some terms will be included here and others 
will not be. You can edit the list of what is checked. Terms retrieving zero 
results are not included in the report.

## Sample report output, HTML and screenshot

Number of studies by audience type. 

Preview the HTML report: http://htmlpreview.github.io/?https://github.com/wendlingd/Browse-InfoSeeking-PubMed/blob/master/InfoSeekingStudies.html

Screenshot: 
<kbd><a href="http://htmlpreview.github.io/?https://github.com/wendlingd/Browse-InfoSeeking-PubMed/blob/master/InfoSeekingStudies.html"><img src="UserStudiesReport.png"></a></kbd>

## Details

Through this basic Python script, the Persons branch of the MeSH tree becomes a
useful tool for accessing research by individual audience types. The
script provides a standing count of studies for each named audience, that you
can retrieve from pubmed.gov.

Run this bibliometric report periodically so you and your staff can have an 
uncomplicated foothold into this type of research.

The script's output is an HTML file with hyperlinks to pubmed.gov, for each Persons 
term that retrieves PubMed records added within the past 6 years.

MeSH changes; you may want to update the csv file here to match the pages
starting from https://www.ncbi.nlm.nih.gov/mesh/68009272. Most yearly updates
are done by January.

## Requirements

- BioPython package, http://biopython.org/DIST/docs/tutorial/Tutorial.html
- MyNCBI account, https://www.ncbi.nlm.nih.gov/books/NBK3842/
