# Fortran programs
This set of Fortran programs were developped to perform or assist with various task related to Single Nucleotide Polymorphism (SNP) data and pedigrees. All are able to handle datasets with several tens of thousands of individuals. 

- `sequoia` : Pedigree reconstruction from SNP data. It has a large overlap with the R package of the same name, for which extensive documentation is available at jiscah.github.io
- `pedigree_checker` : Checks if a pedigree is consistent with the genetic SNP data, by calculating for each offspring - parent - parent trio various probabilities that either or both parents truly are parents, or otherwise related.
- `find_pairs` : Finds in the genetic SNP data all sample pairs which are likely to be duplicates, or parent-offspring, or have a specified other close relationship
- `imputator` : Imputation of missing data, choose from various methods
- `format_geno` : convert genotype files to a different format
- `grm_tool` : Calculate various summary statistics from (potentially huge) GRMs and/or filter out pairs with R values above or below a specified threshold

# Other
Each of the Fortran programs takes genetic SNP data in the same input format as the stand-alone sequoia ( https://github.com/JiscaH/sequoia_notR ). Reformatting from PLINK ped/map or bim/bed/fam formats can be done with the 
format_SNP_data_for_sequoia.sh script. 


## Installation
To run, each program first needs to be compiled by a compiler that supports at least Fortran 2003, e.g. gfortran (https://fortran-lang.org/learn/os_setup/install_gfortran/) or ifort (https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html). Under windows you will also probably need a linux emulator, e.g. cygwin (https://www.cygwin.com/). 

Then, compile with
```
gfortran -O3 <file_name>.f90 -o <program_name>
```

The `imputator` and `sequoia` programs share a set of input/ouput subroutines (`sqa_fileIO.f90`) and general subroutines (`sqa_general.f90`) and can be compiled with

```
gfortran -O3 sqa_fileIO.f90 sqa_general.f90 <file_name>.f90 -o <program_name>
```

where the order of `sqa_fileIO.f90` and `sqa_general.f90` does not matter, but both must come before `imputator.f90` or `sequoia.f90`.  

Program `format_geno.f90` uses only `sqa_fileIO.f90`. 


## Genotype format
The default format used has genotypes coded as 0, 1 or 2 copies of the reference allele, missing values coded as -9, sample IDs in column 1, and no header row. This format can be created as follows:

```
plink --bfile <original_file> --recode A --out TmpFile;
cat TmpFile.raw | tr -s ' ' | cut -d ' ' -f2,7- > sequoia_file.txt;
sed -i '1d;s/NA/-9/g' sequoia_file.txt;
rm TmpFile.*
```
which removes the header, columns 1 (FID) and 3-6 (sex,dam,sire,phenotype) from the .raw file, and replaces all NA's by -9. 

Alternatively, use `format_geno.f90`, which can convert to this format from PLINK .ped/.map files and .raw file, and back again. Programs `imputator` and `sequoia` also support PLINK .ped/.map and .raw files directly. 



## Running
Each program has various command-line options, e.g.

```
./sequoia --geno genotypefile.txt --par --quiet
```

an overview is shown with `./<program_name> --help` . 

For more information, see each program's README file. 

