# Browse InfoSeeking studies in PubMed

## Health-medicine resource for building personas, journey maps, etc.

Browse 'information needs-seeking-use' studies in PubMed with Python-generated 
report

Task: Let's say you are creating a new web resource for a particular 
audience, such as "caregivers." What do you know about the information
needs, information seeking, and information use behaviors of this audience? 
The more you know, the more effective your communication will be. Much as 
been published.

Subject matter experts, product managers, public affairs staff, etc. in
medicine- and health-related disciplines could use an ongoing connection to
this type of research.

Through this basic Python script, the Persons branch of the Medical Subject 
Headings (MeSH) tree, https://www.ncbi.nlm.nih.gov/mesh/68009272, becomes a
useful tool for accessing this research by individual audience types. It
provides a standing count of studies for each named audience, that you can
retrieve from pubmed.gov.

Information introducing the script describes mechanisms, requirements, and 
a few caveats.

Run this bibliometric report periodically so you and your staff can have an 
uncomplicated foothold into this type of research.

## Sample report output

Number of studies by audience type. MeSH page above shows where the indents are.

```
Total records found = 1264
 
To view studies, search at pubmed.gov using "AudienceTermHere"[Mesh] AND 
"Information Seeking Behavior"[Mesh] AND "last 6 year"[dp]
 
Abortion Applicants 0
Adult Children 2
Age Groups 811
Adolescent 286
Adult 752
Child 87
Infant 32
Alcoholics 0
Athletes 4
Caregivers 27
Child of Impaired Parents 2
Child, Abandoned 0
Child, Exceptional 0
Child, Gifted 0
Child, Orphaned 0
Child, Unwanted 0
Consultants 0
Crime Victims 2
Adult Survivors of Child Abuse 0
Criminals 0
Disabled Persons 12
Amputees 0
Disabled Children 4
Mentally Disabled Persons 1
Mentally Ill Persons 0
Persons With Hearing Impairments 2
Visually Impaired Persons 0
Disaster Victims 0
Drug Users 2
Emigrants and Immigrants 12
Undocumented Immigrants 0
Famous Persons 9
Friends 14
Grandparents 0
Homebound Persons 0
Homeless Persons 1
Homeless Youth 0
Jehovah's Witnesses 0
Legal Guardians 1
Proxy 0
Medically Uninsured 4
Men 3
Mentors 2
Minors 0
Missionaries 0
Multiple Birth Offspring 0
Quadruplets 0
Quintuplets 0
Triplets 0
Twins 0
Occupational Groups 200
Administrative Personnel 3
Astronauts 0
Counselors 19
Educational Personnel 4
Emergency Responders 5
Ethicists 0
Farmers 0
Foreign Professional Personnel 2
Government Employees 0
Health Personnel 180
Inventors 0
Laboratory Personnel 0
Lawyers 0
Librarians 3
Military Personnel 5
Miners 0
Pilots 0
Police 3
Religious Personnel 1
Research Personnel 4
Social Workers 1
Parents 90
Fathers 5
Mothers 16
Single Parent 0
Surrogate Mothers 0
Patients 21
Adolescent, Hospitalized 0
Adolescent, Institutionalized 0
Child, Hospitalized 0
Child, Institutionalized 0
Inpatients 3
Outpatients 2
Pedestrians 0
Population Groups 91
Continental Population Groups 59
African Continental Ancestry Group 28
American Native Continental Ancestry Group 7
Asian Continental Ancestry Group 22
European Continental Ancestry Group 19
Oceanic Ancestry Group 3
Ethnic Groups 70
African Americans 23
Amish 0
Arabs 3
Asian Americans 12
Hispanic Americans 27
Mexican Americans 1
Inuits 0
Jews 2
Roma 1
Prisoners 2
Prisoners of War 0
Refugees 2
Research Subjects 1
Healthy Volunteers 0
Sex Workers 1
Sexual Minorities 3
Transgender Persons 1
Sexual Partners 6
Siblings 1
Single Person 0
Slaves 0
Spouses 9
Students 60
Student Dropouts 0
Students, Health Occupations  28
Survivors 29
Adult Survivors of Child Adverse Events 0
HIV Long-Term Survivors 0
Terminally Ill 0
Tissue Donors 6
Blood Donors 2
Living Donors 0
Unrelated Donors 0
Transients and Migrants 4
Transplant Recipients 1
Vegetarians 0
Vegans 0
Veterans 3
Visitors to Patients 0
Volunteers 1
Healthy Volunteers 0
Hospital Volunteers 0
Vulnerable Populations 4
Women 32
Battered Women 0
Dentists, Women 0
Physicians, Women 0
Pregnant Women 28
Women, Working  2
Working Poor 0
```
