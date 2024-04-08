# Manual for `format_geno.f90`


## Description
A Fortran program to convert genotype files to a different format.

See `--help` for command line options


## Supported formats


| Code | Header | ID column # | Columns per SNP | Genotype coding | Missing value |  Notes |
| ---- | ----  | ---- | ---- | ---- |---- | ---- |
| SEQ  | No  | 1 | 1 | 0,1,2 | -9 or -1 | Default format for `sequoia` |
| RAW  | Yes | 2 | 1 | 0,1,2 | NA | PLINK --recode A, columns 1 & 3-6 ignored
| PED  | No | 2 | 2 | 1/2 or A/C/T/G | 0 | paired with .map file, columns 1 & 3-6 ignored 
| LMT  | No | - | 1 | 0,1,2 | not allowed | no space between 'columns', paired with .id file

