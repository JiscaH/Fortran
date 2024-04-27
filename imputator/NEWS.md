## 0.3.7  2024-04-27
- fix --impute-all : imputions for non-genotyped individuals were written to edits
file but not to output genotype. Now added to bottom of genotype file. 
- updated manual


## 0.3.6  2024-04-19
- fix --edits-in: was ignored due to conflict with default options (for do_snpclean)
- fix edits output file: ant_<g> & probs_<g> columns only when method='full';
set prob_<g> to population genotype frequencies for methods het & common


## 0.3.5  2024-04-08
- add 'method' input parameter, with new method 'parent' ('ancestors' = --quick)
- speed increase for methods ancestors and parents
- make sqa_fileIO & sqa_general compatible with both sequoia & imputator


## 0.2.3  2024-02-16
- bugfix in sqa_fileIO.f90 readGeno for inFormat PED
- R script for parallel imputation

## 2024-02-13
- Imputation based on iterative peeling
- Various genotype input & output formats
- Pedigree implemented as a vector of derived type 'individual'


# 0.1.4  2023-12-20
Prototype, imputation based on parent genotypes only