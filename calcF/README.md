# Manual for `calcF.f90`

### A fortran program to calculate pedigree & genomic inbreeding coefficients

Jisca Huisman, 2024-08-11

## Description

This program can

- calculate pedigree inbreeding coefficients, using the method described in Meuwissen & Luo, 1992 (https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-24-4-305)
- calculate genomic inbreeding coefficients
  - correlation between uniting gametes, identical to Fhat3 in GCTA and PLINK (https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM) (`F_uni_u`)
  - the adaptation of the above suggested by Zhang, Goudet, & Weir 2022 (https://www.nature.com/articles/s41437-021-00471-4) (`F_uni_w`) 
- calculate the genomic relatedness between parent pairs
  - can be used to predict the genomic inbreeding coefficients for non-genotyped individuals (e.g. unborn)
  - currently only `R_uni_u` implemented, i.e. identical to that in the GRM returned by GCTA & PLINK
  
  
  
## Requirements & preparation

- a Fortran compiler, e.g. gfortran
  (<https://fortran-lang.org/en/learn/os_setup/install_gfortran/>)

The source code first needs to be compiled with

``` bash
gfortran -O3 sqa_fileIO.f90 calcF.f90 -o calcF
```

where `gfortran` is the name of a Fortran compiler and `-O3` is the
optimisation level. On many platforms `-O4` may be available as well for
a faster program. The choice of program name (after `-o`) is completely
free. `sqa_fileIO.f90` can be found one directory up from here, 
and contains the file input methods shared across several of my other Fortran programs. 


## etc.

To be written. Please run `./calcF --help` to see program options. 