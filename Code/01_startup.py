'''
01_startup.py

https://github.com/wendlingd/Browse-InfoSeeking-PubMed
last updated: 2026-04-04

This script might help sustain a project of this type. You may not need it in your first running
of the primary script.

See requirements.txt for required packages.
'''


#%%
# =========
# Get MeSH
# =========
'''
Collect the MeSH terms that fit your need for information about x group's 
information needs, seeking, and use. The api collection script can use a 
database, or not. For simplicity, this project describes how to use CSV and 
Excel, but the basics of establishing a local SQLite database appear below.

MeSH is a taxonomy of Medical Subject Headings; current information about the 
Persons branch of the MeSH tree is available from the U.S. National Library 
of Medicine: https://meshb.nlm.nih.gov/record/ui?ui=D009272.

I copy what I need and then use the tree notation assign the hierarchy level 
to maintain in the finished HTML report.

- Go to https://meshb.nlm.nih.gov/record/ui?ui=D009272
- Expand every plus sign that you need.
- What I do is copy the whole Person's branch and name it 
locally by the MeSH edition year.
- Or, paste the specific nodes you want, to a text editor.
- In the text editor:
- Replace space-open-square-bracket with “,
- Switch to regular expressions
- Replace \]\r\n with \r\n" (on MacOS you don't need \r)
- Look at result; you may need to Replace \] \r\n with \r\n"
- Replace closing square bracket with nothing
- Add new row 1: 2 column names such as term,tree_number
- Save as CSV

For example, perhaps a person wants:
Occupational groups > Research personnel
Occupational groups > Health personnel (with a few reductions such as only the broad “Nursing” category without specialties)

PubMed searches will 'explode down', meaning, when you use a broader term, narrower 
terms will also be added into the search, so you don't need to include all narrower terms.
'''

import pandas as pd

# Open in Python, whole year or only the terms you need.
broad_mesh = pd.read_csv('2026_mesh_persons.csv')

# Add 'indent' (MeSH hierarchy level) and send back out for manual edits
# The below counts periods, i.e., MeSH tree hierarchy levels. My immediate intention
# here is to plan for HTML indents in the finished report, so users of the HTML page
# will understand categories and sub-categories.
broad_mesh['indent'] = broad_mesh['tree_number'].str.count(r'\.') - 1

# Write out a version you can edit manually, that later will become your API source dataframe
broad_mesh.to_csv('usable_mesh.csv', index=False)

'''
Name your result Data/Processed/current_project_mesh.xlsx (for example)
Background: Excel can be useful for reducing the rows needed and then 
running them against PubMed to experiment with results. While selecting,
I compare back and forth with https://meshb.nlm.nih.gov/record/ui?ui=D009272.
'''


#%%
# =====================================
# Optional: Customize the search query
# =====================================
'''
Alter the base search strategy to suit your need. A resource that may or may not facilitate
your work with the PubMed interface: URL encoder/decoder at 
https://meyerweb.com/eric/tools/dencoder/
My default search strategy has been:

(NAMED AUDIENCE)[Mesh]
AND ("Information Seeking Behavior"[Mesh] OR "user stud*"[Title/Abstract] OR "user research"[Title/Abstract])
AND ("last 6 years"[PDat])

The following were tested as possible additions to the strategy, but not all results
appeared relevant for the current purpose. Consider adding if appropriate for your task.

- usability[Title/Abstract]
- "User-Centered Design"[Mesh]
- "Human-Centered Design"[Title/Abstract]
- "Consumer Behavior"[Mesh]
- "Data Collection"[Mesh]
- "Qualitative Research"[Mesh]
- "case stud*"[Title/Abstract]
- "interview*"[Title/Abstract]
- "Information Literacy"[Mesh]
'''


#%%
# ================================
# Optional: create local database
# ================================
'''
Not required. Potentially useful for long-term projects.
'''

import os
from pathlib import Path
import sqlite3

from pathlib import Path
# Build path in a cross-platform way
target_dir = Path.home() / "Documents" / "Browse-Infoseeking"
# Change working directory
os.chdir(target_dir)


# Create the Data/Local directory if it does not already exist
db_dir = Path("Data") / "Local"
db_dir.mkdir(parents=True, exist_ok=True)

# Database file
db_path = db_dir / "mesh_counts.db"


with sqlite3.connect(db_path) as conn:
    cur = conn.cursor()

    # Create table to store all MeSH terms
    cur.execute("""
    CREATE TABLE IF NOT EXISTS persons2026 (
        mesh_id INTEGER PRIMARY KEY AUTOINCREMENT,
        term TEXT NOT NULL,
        count INTEGER,
        tree_number TEXT
    )
    """)

    # Set AUTOINCREMENT starting point so first generated mesh_id is 100
    cur.execute("""
    INSERT OR REPLACE INTO sqlite_sequence (name, seq)
    VALUES ('persons2026', 99)
    """)


# Create table to store API results
from pathlib import Path
import sqlite3

db_path = Path("Data") / "Local" / "mesh_counts.db"
db_path.parent.mkdir(parents=True, exist_ok=True)

with sqlite3.connect(db_path) as connection:
    cursor = connection.cursor()

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS api_result (
        api_result_id INTEGER PRIMARY KEY AUTOINCREMENT,
        indent INTEGER,
        mesh TEXT,
        query TEXT,
        pubmed_count INTEGER,
        tree_number TEXT,
        status TEXT,
        date_run TEXT
    )
    """)

    # Set the next generated ID to 100 if the table is brand new
    cursor.execute("""
    INSERT OR IGNORE INTO sqlite_sequence (name, seq)
    VALUES ('api_result', 99)
    """)


print(f"Database created/updated at: {db_path}")
print("Tables `persons2026` and `api_result` are ready for use.")


#%%
