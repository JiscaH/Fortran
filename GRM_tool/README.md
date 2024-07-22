# Manual for `grm_tool.f90`

### A fortran program to summarise, subset, and filter GRMs

Jisca Huisman, 2024-07-22

## Description

This program can

- calculate various summary statistics from (potentially huge) Genomic
  Relatedness Matrices (GRMs);
- subset the GRM based on a list with (pairs of) individuals;
- filter out pairs with relatedness (R) values above, below or between
  specified thresholds;
- filter out individuals with inbreeding values (GRM diagonal) above,
  below or between specified thresholds

Additional functionality may be added in the future, please feel free to
contact me with suggestions.

Computational time may be more than an hour for very large GRMs, which
is predominantly spent on unzipping and loading the data.

## Preparations

### Requirements

- a Fortran compiler, e.g. gfortran
- gzip (or pigz), to decompress the .grm.gz file containing the GRM

### Compilation

Run `gfortran -O3 grm_tool.f90 -o grm_tool` (where you may substitute
gfortran by another compiler of your choice).

## Input files

- `--in` : input file pair with GRM, as returned by
  `plink --make-grm-gz`
  (<https://www.cog-genomics.org/plink/1.9/distance>) or
  `gcta --make-grm`
  (<https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM>).
  Extensions .grm.id and
  .grm.gz[^](except%20when%20option%20%60--notgz%60%20is%20used,%20then%20extension%20.grm%20is%20added%20and%20the%20file%20is%20presumed%20to%20be%20a%20plain%20text%20file)
  are added. Binary files (.grm.bin) are not currently supported.
- `--only` : individual subset, only pairs where either individual is
  listed are considered. A text file with no header, and either two
  columns with IDs in the second column (same format as for
  `plink --keep`), or a single column with IDs (any additional columns
  are ignored if there are \>2 columns).
- `-only-among` : as `--only`, but only pairs where *both* individuals
  are listed are considered.
- `--only-pairs`: A text file with no header and 2 columns: ID1 and ID2
  of each pair.

## Program options

### Subset

There are 3 ways of subsetting the GRM using lists of IDs: -
`--only-pairs`: only the specified pairs. The order of the pair (A-B vs
B-A) is automatically corrected to match the .grm.gz file, which by
default only contains the lower-triangular GRM - `--only-among`: only
pairs for which BOTH individuals are listed in the file with IDs -
`--only`: only pairs for which EITHER individual is listed in the file
with IDs

### Summary statistics

Currently, the only command line options for the summary statistics are
to turn calculation off (`--no-summary`) and to specify the outfile
(`--summary-out`). Turning the calculation off potentially allows
filtering of GRMs with more individuals, see section ‘Method’.

The summary output currently looks as follows:

              group           part             n_pairs         n_snps     minimum        maximum        mean           std_dev           count_<-0.5         count_<0.625        count_>0.875        count_>1.25 
             total        diagonal                 200         495.07       0.000000       1.265416       0.984194       0.077719                   0                   0                 189                   1
             total        between                19900         490.21      -0.221723       0.756767      -0.004994       0.101084                   0               19892                   0                   0
             across       diagonal                   4         496.50       0.000000       1.087468       1.040321       0.046471                   0                   0                   4                   0
             across       between                  790         491.60      -0.178162       0.567930      -0.005751       0.069678                   0                 790                   0                   0
             among        diagonal                   4         496.50       0.000000       1.087468       1.040321       0.046471                   0                   0                   4                   0
             among        between                    6         493.00      -0.006430       0.138561       0.055076       0.053186                   0                   6                   0                   0
             pairs        diagonal                   1         498.00       0.000000       0.977654       0.977654            NaN                   0                   0                   1                   0
             pairs        between                    3         491.33      -0.006430       0.098249       0.044369       0.052408                   0                   3                   0                   0

where the groups are

- *total*: in the full GRM
- *across*: between members of the `--only(-..)` list and all other
  individuals, AND among members of the list
- *among*: among members of the `--only(-..)` list; for `--only-pairs`
  this is among all possible pairwise combinations of the individuals
  included in the pairs
- *pairs*: Only for `--only-pairs`: just among the specified pairs.

and the ‘part’ column separates diagonal elements of the GRM (= 1+F) and
the between-individual off-diagonal elements

When the program is run without any `--only(-..)`, just the ‘total’
group is returned; the ‘pairs’ rows are only returned for `--only-pairs`

The count thresholds are currently hardcoded. To change, search for
‘sumstats_counts’ in the source code; do not forget to also change the
output column header accordingly, via ‘sumstat_lbls’. Please feel free
to add more summary statistics, or contact me for assistance if you are
not familiar (enough) with Fortran to make the necessary edits yourself.

Note that when reading this file into R, its default behaviour is to
change the ‘\<’, ‘\>’ and ‘-’ in the column names into a dot. To prevent
this, use `read.table(file, header=TRUE, check.names=FALSE)`.

### Filtering

A subset of the GRM can be written to a text file based on their R
values, e.g. to check any weird values (e.g. R \> 1.5 or R \< -1.0) or
to find off-diagonal pairs which may be 1st degree relatives but not
duplicates or inbred (e.g. 0.35 \< R \< 0.65).

> [!Note] 
> This is meant as a tool for data quality control, and I strongly advise against using this as primary method to identify 
> relatives. For that, please use sequoia, `find_pairs.f90` (both in
> <https://github.com/JiscaH/Fortran>), or other likelihood-based
> methods.

Separately for the diagonal and off-diagonal (‘between’) elements, one
or two thresholds can be set:

- `diag-lower`, `betw-lower` : the lower threshold, pairs with values
  above this thresholds are written to the output file;
- `diag-upper`, `betw-upper` : the upper threshold.

One threshold may be left unspecified, and the `upper` threshold may be
higher or lower than the `lower` threshold, giving four different
possibilities as illustrated as A–D in the figure (grey = selected pairs
written to the output).

<img src="grm_tool_manual_files/figure-gfm/filter_illustration-1.png" width="60%" />

If combined with `--only`, the subset is within the specified group
only. So for example,

- `--only those.txt --betw-lower 0.4` returns pairs with R\>0.4 for
  which one or both individuals are listed in ‘those.txt’
- `--only-among those.txt --diag-upper 0.9` returns individuals in
  ‘those.txt’ with GRM diagonal \< 0.9
- `--only-pairs those_pairs.txt --betw-upper 0.2` returns pairs in
  ‘those_pairs.txt’ with R\<0.2

#### Infinity

> [!Note] 
> Filtering will currently ignore any values of `+Inf` or
> `-Inf`, as the thresholds do not default to infinity but to
> `HUGE(0D0)`, which is the largest ‘double precision’ value that
> Fortran can store. Whether those are present in the data can be
> checked in the summary statistics as

``` fortran
write(42,*)  'count_inf ',  COUNT(GRM >= HUGE(0D0))
```

If this is an issue, please let me know as it would be fairly
straightforward to change.

### Histogram

Counts per histogram bin are calculated with `--hist`, with separate
counts for the diagonal and off-diagonal (‘between’) elements. If
combined with `--only(-..)`, the groups are identical to those for the
summary: ‘total’, ‘across’, ‘among’, and ‘pairs’.

The result is written to a 4-column text file (default hist_counts.txt,
specify with `--hist-out`), with columns:

- group: as for the summary: ‘total’, ‘across’, ‘among’, and ‘pairs’
- part: diagonal vs between (off-diagonal)
- lower_bound: lower bound of each bin.
- count: number of pairs with R \> \[this lower_bound\] and R \<= \[next
  lower bound\]

> [!Note] 
This file format differs from the January 2024 program
version.

The lower bounds are equidistant and default from -1.5 to +2.0 with a
stepsize of 0.05. These can be adjusted as optional arguments to
`--hist`, in which case all 3 must be provided, in the order first,
last, step.

The output is (should be) identical to the output of the R command

``` r
  table(cut(GRM, breaks=seq(first, last, step)))
```

All bins are closed on the right and open on the left, any values \<=
first or \> last are discarded, with a warning. It might therefore be
advisable to edit the summary statistics to count the number of values
on and off the diagonal exceeding ‘last’. Note that the warning will be
given only after all the data has been processed.

The output can e.g. be visualised in R as follows:

``` r
H <- read.table('my_grm_hist_counts.txt', header=TRUE)
library(magrittr)
H %>% dplyr::filter(group=='total' & part=='diagonal') %>%
  plot(count ~ lower_bound, data=., type='h', xlim=c(-.4, 1.5),
       main='Diagonal', ylab='Count', xlab='R' )
# - or - 
H %>% dplyr::filter(group=='total' & part=='diagonal') %>%
  barplot(count ~ lower_bound, data=., space=0, 
          main='Off-diagonal', ylab='Count', xlab='R')
# - or -
HH_list <- plyr::dlply(H, c('group', 'part'), function(Hs) {
  tmp <- list(breaks = c(Hs$lower_bound,
                         Hs$lower_bound[nrow(Hs)] +
                           (Hs$lower_bound[2]-Hs$lower_bound[1])),
              counts = Hs$count,
              xname = 'R')
  class(tmp) <- 'histogram'
  return(tmp) })
plot(HH_list[['among.between']], xlim=c(0, 1.3))
```

## Method

The compressed .grm.gz file is decompressed and processed as a
continuous stream of data via a so-called ‘named pipe’, as described at
<https://genomeek.wordpress.com/2012/05/07/tips-reading-compressed-file-with-fortran-and-named-pipe/>
. By default `pigz` is used for decompression, but this can be changed
to gzip with `--zipper gzip`. For other programs, change the following
line in the source code:

``` fortran
call EXECUTE_COMMAND_LINE("(pigz -dc  "//trim(grmFile)//".grm.gz > grmpipe) &")
```

The speed of the program is limited by the decompression speed. Whereas
compression of data can be done in parallel on multiple threads,
decompression cannot (according to <https://zlib.net/pigz/pigz.pdf> ).

For the filtering, the data is taken from the stream, checked against
the specified thresholds, written to file (or not), and discarded. The
amount of computer memory required by the program does thus not increase
with the size of the GRM, although the output file may get very large
depending on the thresholds used.

For the summary statistics and `--hist`, after 5%, 10%, …, of the data
is read in (default `--chunks` is 20), summary statistics are calculated
and memory space reused for the next 5% of data. After all data are
processed, the sum/min/max/mean is calculated over the 21 data chunks
(20 equal sized + left overs). For the standard deviation, the
‘extension to K groups’ at
<https://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation?noredirect=1&lq=1>
is used.

When the number of individuals is large relative to the amount of
available memory, a larger number of chunks can be specified with with
`--chunks <x>`. Limited testing suggests computational time is minimised
at an intermediate number of chunks.

## Example data

A small mock GRM is available to try out this program and its various
setting. The files ‘griffin.grm.gz’ and ‘griffin.grm.id’ contain a GRM
generated from SNP data simulated from ‘Ped_griffin’ in the sequoia R
package (<https://CRAN.R-project.org/package=sequoia>). It includes 200
individuals, including 2 inbred pairs with R\>0.75.
