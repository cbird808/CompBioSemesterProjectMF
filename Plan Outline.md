Plan outline

Write a bash script that extracts the file name, read number, and first 16 base pairs of “.fq” files for both treatments 
(albatross and contemporary) within a directory and create a “.csv” file which can be fed into R. Within the script, 
R will be used to determine whether each base pair matches the expectation for that location, create summary statistics, 
and generate graphs representing the percentage of correct matches for each base pair separate for each barcode and treatment. 
The goal is to see whether there is a difference between the star activity present in albatross vs contemporary samples.

Timeline
1.	(complete) Create script to extract filenames read number, and first 16 base pairs of “.fq” files for both treatments (albatross and contemporary) within a directory and create a “.csv” file which can be fed into R. 
2.	(by Oct 11th) Write R script to determine whether each base pair matches the expectation for that location, create summary statistics, and generate graphs representing the percentage of correct matches for each base pair separate for each barcode and treatment.
3.	(by Oct 31st) Create shell script to contain both the bash script and the R script necessary to make this function at the command line.
4.	(by Nov 30th) Finalize script and make report in markdown
