'''
02_InfoSeeking-studies-by-audience.py

https://github.com/wendlingd/Browse-InfoSeeking-PubMed
last updated: 2026-04-04

This script uses the eutils API to pull descriptor counts from PubMed.
It's Python for browsing the MeSH 'Persons' branch in PubMed.

Follow the documentation:
 - MyNCBI: https://www.ncbi.nlm.nih.gov/myncbi/
 - E-utilities Quick Start: https://www.ncbi.nlm.nih.gov/books/NBK25500/
'''


#%%
# ========
# Startup
# ========

import os
import json
from Bio import Entrez
import pandas as pd
import time
# import requests
import urllib.parse
import sqlite3 # optional

import datetime
from datetime import date, datetime, timedelta

from pandas import read_excel

# from datetime import timedelta
today = datetime.today().strftime('%Y-%m-%d')

from pathlib import Path
# Build path in a cross-platform way
target_dir = Path.home() / "Documents" / "Browse-Infoseeking"
# Change working directory
os.chdir(target_dir)

# Current working dir
print(f'Current working directory: {os.getcwd()}')


#%%
# =======================
# Set variables and test
# =======================
'''
Derive a start date. Maybe the past 6 years for your first run? And later, 
expand and save locally; go backward as you accumulate audience knowledge.
'''

# Load NCBI credentials
with open("../creds/pm.json", "r", encoding="utf-8") as f:
    creds = json.load(f)

Entrez.email = creds["email"]
Entrez.api_key = creds["key"]
# Optional but recommended by NCBI / Biopython
Entrez.tool = "pubmed_count_collector app"

# Test your authentication by asking for the list of supported databases
with Entrez.einfo() as handle:
    record = Entrez.read(handle)

print(f'Can you authenticate? Test result: \n{record["DbList"]}')


#%%
# =======================
# Bring in audience list
# =======================
'''
Usually MeSH is updated near the beginning or end of the calendar year.
'''

# Load the user groups file - selected levels and terms from the Persons branch (M01)
# of NLM's MeSH vocabulary.
# audienceList = pd.read_csv('Data/MatchFiles/current_project_mesh.csv')
audience_list = pd.read_excel("Data/Processed/current_project_mesh.xlsx")

# Limit when testing
# audience_list = audience_list.sample(5)

print(f'audienceList shape is {audience_list.shape}\n')


#%%
# =============================
# Collect data into api_result
# =============================
'''
Interact with API (Entrez.esearch).
'''

results = []
row_count = str(len(audience_list))

print(f'\nProcessing {row_count} items...')

for i, (_, row) in enumerate(audience_list.iterrows(), start=1):
    current_term = row["term"]
    indent = row["indent"]
    tree_number = row["tree_number"]

    # Skip blanks / nulls
    if pd.isna(current_term) or str(current_term).strip() == "":
        results.append({
            "indent": indent,
            "mesh": current_term,
            "query": None,
            "pubmed_count": None,
            "tree_number": tree_number,
            "status": "blank MeSH"
        })
        continue

    mesh_term = str(current_term).strip()

    search_strategy = f'''("{current_term}"[mh]) AND ("Information Seeking Behavior"[mh] OR "information seeking"[ti] OR "information needs"[ti] OR "user stud*"[tiab] OR "user research"[tiab]) AND ("last 6 years"[PDat])'''.strip()

    try:
        with Entrez.esearch(
            db="pubmed",
            term=search_strategy,
            retmode="xml",
            retmax=0
        ) as handle:
            result = Entrez.read(handle)

        count = int(result["Count"])

        results.append({
            "indent": indent,
            "mesh": current_term,
            "query": search_strategy,
            "pubmed_count": count,
            "tree_number": tree_number,
            "status": "ok"
        })
        print(f"{i} of {row_count}: {current_term}")

        time.sleep(1)

    except Exception as e:
        results.append({
            "indent": indent,
            "mesh": current_term,
            "query": search_strategy,
            "pubmed_count": None,
            "tree_number": tree_number,
            "status": f"error: {e}"
        })

api_result = pd.DataFrame(results)

# ----------------------
# Light post-processing
# ----------------------
# Add today's date
api_result['date_run'] = today

# Ensure MeSH Tree order
api_result = api_result.sort_values(by=['date_run', 'tree_number'], ascending=[False, True])

# ----------------------
# api_result vs. report
# ----------------------
'''
Use api_result if you want to analyze what retrieved zero records; perhaps you will want
to change the search strategy, adding terms or changing the time horizon. api_result 
will be written to the database (if a db used).
'''

# For the HTML report, drop rows with no PubMed records
report = api_result.loc[api_result['pubmed_count'] != 0].copy()

print(f'\napi_result is {api_result.shape}')
print(f'report is {report.shape}')


#%%
# ==========================
# Optional: Write to SQLite
# ==========================
'''
See 01_startup.py for info. In mine I include zero counts.

Default update is to replace previous content.
'''

db_path = Path("Data") / "Local" / "mesh_counts.db"

# Keep only the columns that belong in the table, excluding api_result_id
api_result_to_load = api_result[
    ["indent", "mesh", "query", "pubmed_count", "tree_number", "status", "date_run"]
].copy()

with sqlite3.connect(db_path) as connection:
    api_result_to_load.to_sql(
        "api_result",
        connection,
        if_exists="append",
        index=False
    )

print(f"{len(api_result_to_load):,} rows replaced in `api_result`.")


#%%
# ===============
# Output to HTML
# ===============
'''
Put the data into an easy-to-use, HTML report format.
'''

# Limit when testing
# report = report.sample(5)

# Keep only the columns you want, in the order you want
report = report[['indent', 'mesh', 'query', 'pubmed_count']].copy()

def make_link(query: str, count: int) -> str:
    """Convert a PubMed query and count into an HTML hyperlink."""
    encoded_query = urllib.parse.quote_plus(query)
    url = f"https://pubmed.ncbi.nlm.nih.gov/?term={encoded_query}"
    return f'<a href="{url}" target="_blank">{count:,}</a>'

def make_html_row(row):
    count = int(row["pubmed_count"]) if row["pubmed_count"] is not None else 0
    count_link = make_link(row["query"], count)
    return f'<li class="indent{row["indent"]}">{row["mesh"]} ({count_link})</li>'

# Generate the <li> rows
report["HTML"] = report.apply(make_html_row, axis=1)

# HTML skeleton
html_builder = [
    "<html><head><title>Audience-specific information studies at pubmed.gov</title>",
    "<style>",
    "  body {margin-left:1em; font-family: Arial, Helvetica, sans-serif;}",
    "  ul {list-style-type: none; margin-left:0; padding-left:0;}",
    "  .indent1 {padding-left:0em;}",
    "  .indent2 {padding-left:2em;}",
    "  .indent3 {padding-left:4em;}",
    "  .indent4 {padding-left:6em;}",
    "  .indent5 {padding-left:8em;}",
    "</style></head><body>",
    "<h1>Understand your audiences: Health-medical information-seeking-behavior studies at pubmed.gov</h1>",
    f"<p><strong>{today}</strong> &ndash; Use this report to understand information-needs, information-seeking, and ",
    "information-use behaviors for various audiences for biomedical information. Example: A link for &quot;caregivers&quot; ",
    "would run the following search strategy at ",
    '<a href="https://pubmed.ncbi.nlm.nih.gov">https://pubmed.ncbi.nlm.nih.gov</a>: Caregivers[Mesh] AND ',
    '("Information Seeking Behavior"[Mesh] OR "information seeking"[Title] OR "information needs"[Title] OR '
    '"user stud*"[Title/Abstract] OR "user research"[Title/Abstract]) AND "last 6 years"[PDat]. This report ',
    "uses a <strong>partial and selected list</strong> of terms from the ",
    '<a href="https://meshb.nlm.nih.gov/record/ui?ui=D009272">Persons branch of the MeSH tree.</a> Terms with zero ',
    "results are not shown. Counts will go out of date. ",
    '<a href="https://github.com/wendlingd/Browse-InfoSeeking-PubMed">https://github.com/wendlingd/Browse-InfoSeeking-PubMed.</a></p>',
    "<ul>"
]

# Append all list items
html_builder.extend(report["HTML"].tolist())

# Close the HTML
html_builder.extend([
    "</ul>",
    "</body>",
    "</html>"
])

# Write the output file
output_path = "Docs/Reports/InfoSeekingStudies.html"
with open(output_path, "w", encoding="utf-8") as f:
    f.write("\n".join(html_builder))

print("HTML report written to Docs/Reports/InfoSeekingStudies.html")


#%%

