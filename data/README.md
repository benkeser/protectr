# Data specification

This README outlines formatting required for compatibility with the code in this repository. 

In the source data, each client should have one row of data per encounter.

There are two types of variables -- those that vary by encounter and those that are constant across visits.

## Static variables across visits

These variables should have the same value for each row of data for a client.

- `id` = numeric variable identifying each client
- `enroll_date` = date variable identifying date of enrollment into care
	- confirm timeline for enrollment into care (e.g., > 01/01/2018?)
	- should not be `NA` for anyone
- `art_start_date` = date of ART initiation
	- should it always be > `enroll_date`?
	- should be `NA` if client not observed to start ART
- `tpt_start_date` = date of first TPT initiation after enrollment into care
	- should always be > `enroll_date`
	- should be `NA` if client not observed to start ART
- `tb_diagnosis_date` = date of first TB diagnosis after enrollment into care
	- should always be > `enroll_date`
	- should be `NA` if client not observed to have TB diagnosis after enrollment into care
	- ultimately, we exclude individuals with TB within 30 days of enrollment into care
		- this should happen upstream for efficiency of sensitivity analyses
- `right_cens_date_tb` = date variable identifying the last date at which the TB status of the client could be ascertained with reasonable certainty
	- for individuals who are never lost to follow-up for a period longer than 6 months, this should be the date of the last clinic visit
	- for individuals who are lost to follow-up for a period longer than 6 months, this should be the date of the last clinic visit prior to becoming lost to follow-up
	- should not be `NA` for any client
- `last_visit_date` = last date that client was seen in the clinic
	- should not be `NA` for any client
- `cd4_count_date` = date on which initial CD4 count was measured
	- can be `NA` if no CD4 count was taken
- `cd4_count` = value for initial CD4 count 
	- can be `NA` if no CD4 count was taken
- other baseline variables (e.g., sex, age) 
	- should be in form appropriate for modeling using glm

## Visit-level variables

- `visit_date` = date variable identifying the date of each encounter included in the data
	- the first `visit_date` should = `enroll_date`
- `art_adherence` = a numeric variable indicating level of ART adherence
	- may need to assign arbitrary numeric values to categorical labels
	- can be missing for some visits
	- should be `NA` if individual has not yet initiated ART 
- `tb_sx` = a numeric variable indicating TB symptom screen results at each encounter
	- may need to assign arbitrary numeric values to categorical labels
	- can be missing for some visits
