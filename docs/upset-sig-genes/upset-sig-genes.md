## upset-sig-genes.R

This script produces a basic UpSet plot from a list of DESeq2 sig files.
For example:

``` bash
cd test_data
Rscript ../upset-sig-genes.R \
--plot_file upset-plot-test.pdf \
--output_file upset-plot-test-out.tsv \
set1.sig.tsv set2.sig.tsv set3.sig.tsv
```

![UpSet plot of 3 sets, showing the sizes of each set and the sizes of
the intersections between the
sets](upset-plot-test.png "UpSet plot of 3 sets")

The script tries to automatically label the set using the file names. If
the files are name ‘sig.tsv’ and in separate directories, the directory
name is used (after removing deseq2- from the beginning if it is there).
Otherwise if the files are named \*.sig.tsv, then .sig.tsv is removed
from the end and the remainder is used as the set labels.

For example, inputs files named

deseq2-hom-vs-wt/sig.tsv  
deseq2-hom-vs-het/sig.tsv  
deseq2-het-vs-wt/sig.tsv

will be labelled hom-vs-wt, hom-vs-het and het-vs-wt

And files named

deseq2/hom-vs-wt.sig.tsv  
deseq2/hom-vs-het.sig.tsv  
deseq2/het-vs-wt.sig.tsv

will also be labelled hom-vs-wt, hom-vs-het and het-vs-wt

It also outputs a file of which genes are in which intersections. The
intersections are numbered by converting the binary representation of
which sets are included into a decimal. For example, if there are 3 sets
(set1, set2 and set3), the intersection numbered 5 represents the
elements in both set 1 and set3, but not set2

    set1 set2 set3
     1    0    1   = 5

The sets are order by size, which means they may appear in a different
order to the command and this includes the naming of the intersections.
If you wish to keep the sets in the same order as that of the command
set the `--keep_order` option. In this case, the sets will also be in
this order in the UpSet plot, with the first set at the bottom of the
y-axis.

``` bash
Rscript ../upset-sig-genes.R \
--plot_file upset-plot-test-ordered.pdf \
--output_file upset-plot-test-ordered-out.tsv \
--keep_order \
set2.sig.tsv set3.sig.tsv set1.sig.tsv
```

![UpSet plot of 3 sets, showing the sizes of each set and the sizes of
the intersections between the
sets](upset-plot-test-ordered.png "UpSet plot of 3 sets not ordered by set size")
