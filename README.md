# metaDRMetRap

Data metadata for "Rapamycin, not metformin, mirrors dietary restriction-driven lifespan extension in vertebrates: a meta-analysis"
Paper here: https://onlinelibrary.wiley.com/doi/10.1111/acel.70131  
Data and Code: https://zenodo.org/records/15673918  

Edward R. Ivimey-Cook*@, Zahida Sultanova*@, and Alexei A. Maklakov  
*these authors contributed equally  
@corresponding authors: e.ivimeycook@gmail.com and zahida.sultanova@uea.ac.uk  

Data: analysis_data.csv  
Code: create Code folder and run script "Analysis_Script.R". All other scripts are sourced within this one.  

Data columns (all names are cleaned when loading):  
id = 1:911 unique individual ID  
author list: list of authors  
title: title of paper  
Year: year paper was published (extracted from Scopus/WoS)  
m_control: control lifespan measure  
n_control: sample size for control  
sd_control: standard deviation for control (left blank if not available)  
se_control: standard error for control (left blank if not available)  
m_treatment: treatment lifespan measure  
n_ treatment: sample size for treatment  
sd_ treatment: standard deviation for treatment (left blank if not available)  
se_ treatment: standard error for treatment (left blank if not available)  
m_Measure: median or mean lifespan  
m_raw: Whether raw data was given to calculate measures.  
m_time: unit of lifespan measure  
m_Location: Location of lifespan measure in the paper.  
Species: species sampled  
Treatment_Type: Type of treatment used (e.g. fasting, percent reduction, rapamycin..)  
Treatment_Overall: Broad treatment groups (Rapamycin, Metformin, DR).  
Treatment_Name: Treatment name in paper  
Sex: Sex studied  
Strain: The name of any strain used if appropriate.  
Other Variables: Whether any other variables tested in the experiment.  
Notes: Additional notes  
mice900_keep: whether or not the mice paper passed the 900 day rule.  

R version and packages loaded (with versions): see session.txt file.   




