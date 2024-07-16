## General
This Fortran program performs:

- `--dup`: identification of pairs of likely duplicated samples (or monozygotic twins)
- `--par`: assignment of genotyped parents to genotyped offspring
- `--ped`: full pedigree reconstruction, including clustering of siblings sharing 1 or 2 non-genotyped parents, and grandparent assignment
- `--maybePO`, `--maybeRel`: find potential (remaining) parent-offspring or related pais
- `--pairs`: calculate for the specified pairs the likelihoods for 7 different relationships

This Fortran program is equivalent to the R package https://CRAN.R-project.org/package=sequoia, except that the input & output does not go via R but via text files and command-line arguments. 
Very large datasets can be incredibly slow or even impossible to load into R; this non-R version works for up to 50K - 100K genotyped individuals. 

Please see the manual at
https://github.com/JiscaH/Fortran/blob/main/sequoia/Sequoia_stand-alone_manual.pdf



## Compilation
A Fortran compiler is needed, such as e.g. gfortran. In Windows you also need a linux emulator, e.g. cygwin (although other approaches may be possible).
```
gfortran -O3 Sequoia_SA.f90 -o sequoia
```

> [!NOTE]
> previous versions relied on `-std=95 -fall-intrinsics`, these flags are no longer necessary and will in fact cause compilation to fail)


## File formats

Template files can be generated using the R package `sequoia`:
```
library(sequoia)
writeSeq(SeqList = SeqOUT_griffin, GenoM = Geno_griffin, 
         OutFormat = "txt", folder = "SequoiaFileTemplates")
```

### Genotype file

Four different genotype formats are supported:

- `--genoformat PED`: PLINK .ped/.map file pair https://www.cog-genomics.org/plink/1.9/formats#ped 
- `--genoformat RAW`: from PLINK `--recode A`, .raw  https://www.cog-genomics.org/plink/1.9/formats#raw
- `--genoformat SEQ` (default): as .raw with SNP genotypes coded as 0/1/2, but no header, IDs in column 1 and no additional non-genotype columns, and missing values coded as a negative integer rather than `NA`. 
- `--genoformat LMT`: .geno/.id file pair; No header row, genotypes coded as 0/1/2 without spacing, missing values not allowed, IDs in separate file.

The default `SEQ` format is fastest to read, as it does not require any internal conversion. It can be created e.g. from standard PLINK bed/bim/fam format as follows:
```
plink --bfile FileNameIN --recode A --out FileNameOUT
mv FileNameOUT.raw tmp.raw
cat tmp.raw | tr -s ' ' | cut -d ' ' -f2,7- > FileNameOUT.raw
sed -i '1d' FileNameOUT.raw 
sed -i 's/NA/-9/g' FileNameOUT.raw
rm tmp.raw
```
which

- recodes to 1 column per SNP, coded as 0/1/2 copies of minor allele, missing = NA
- drops columns 2--6 (family ID, sex, phenotype, etc)
- drops header row
- replaces missing value code NA by -9

### Lifehistory file

A file with

- ID; order may differ from genotype file
- Sex; 1=female, 2=male, 3=unknown, 4=hermaphrodite
- birth 'year': discrete time unit of birth; parents are never born in the same time unit as their offspring. Missing = any negative NUMERIC value
- earliest possible birth year, if exact unknown (optional column)
- latest possible birth year, if exact unknown (optional column)
- latest possible year of reproduction (time unit in which offspring is born; can differ from year of death e.g. for males in mammals or in plants)


### Ageprior file

In its simplest form specifies the minimum and maximum age of parents, but can also be used to convey more detailed information about the age distribution of parents and between siblings. See https://jiscah.github.io/articles/vignette_age/book/index.html


### Additional files
Please see the PDF manual.


## Running
Parameter values are stored in `SequoiaSpecs.txt`. Many of these can be overruled on the command line, as described in the PDF manual. 
```
./sequoia --geno geno.txt --par --ped --verbose
```

All command line options are shown with
```
./sequoia --help
```


## Version
Its version will typically be functionally identical to the latest R package beta version. To obtain an older version, please send an email. 
