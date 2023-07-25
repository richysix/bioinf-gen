# Plot counts from transcriptomic data

[Go to script](../../graph_counts_by_group_facet.R)

This is a script to produce counts plots from RNA-seq (or DETCT) data.  
For information on all the options run

    Rscript graph_counts_by_group_facet.R --help

The required arguments are a samples file and a sig file. It expects the
sample file to have a header with the column names. It expects one of
the columns to be called “sample”. e.g.

    sample  condition
    sample_1    wt
    sample_2    mut

There is an example samples file and sig file in the test_data directory
of this repository. The simplest way to run the script would be this:

``` bash
Rscript graph_counts_by_group_facet.R \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv
```

![Basic count plot showing the normalised counts for the wt condition as
blue circles and the mutant condition as orange
circles](test_counts_basic.1.png "Basic count plot")

By default it tries to use a column named condition in the samples file
as the x variable. The default is to colour the points by condition as
well.

The --x_variable option allows you to name a column to use as the x axis
and you can also specify columns to use for colours (--colour_variable)
and shapes (--shape_variable). The --crossbar option plots a bar to show
the mean/median of each group. It is also possible to supply a file of
Ensembl gene ids and the script will only make plots for those gene ids
in the data.

The example below uses almost all the available options

``` bash
# create a test ids file
echo "ENSTEST006" > test_gene.txt

graph_counts_by_group_facet.R \
--output_file test-condition-sex-treatment.pdf \
--genes_file test_gene.txt \
--x_variable condition \
--colour_variable condition \
--colour_palette 'wt=#0000ff,mut=#ff0000' \
--shape_variable sex \
--facet_variable treatment \
--crossbar median \
--width 12 \
--height 8 \
--theme_base_size 14 \
--rotate_xaxis_labels \
--seed 7635 \
--no_pvalue \
test_samples.tsv test_rnaseq_data.tsv
```

![Count plot showing the normalised counts for the wt condition in blue
and the mutant condition in red. The points are split by the control or
treated. Sex is displayed as circle for female and squares for
male](count_plot_everything.1.png "Count plot by condition by treatment")

There is also an option to supply a file of p-values which be added to
the pot with a line from one condition to another. The file should look
like this:

| GeneID     | condition1 | condition2 |    adjp | log2fc |
|:-----------|:-----------|:-----------|--------:|-------:|
| ENSTEST001 | wt         | mut        | 6.52e-5 |  3.378 |

The `log2fc` column is optional. The --asterisks option converts the
p-values to asterisks.

``` bash
Rscript graph_counts_by_group_facet.R \
--pvalue_file rnaseq_adjusted_pvalues.tsv --asterisks \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv
```

![Count plot showing the normalised counts for the wt condition as blue
circles and the mutant condition as orange circles. A line with three
asterisks above it, shows that the difference between the two groups is
significant.](test_counts_pval.1.png "Count plot by condition with p-values")

The other options are:

- --log10 - Use a log10 scaled y-axis
- --output_data_file Output an Rdata file of the plot objects
- --no_jitter - removes the jitter from the points
- --seed - random seed to make the jitter reproducible
- --no_pvalue Don’t add a pvalue to the plot title
- --detct - input data is DeTCT rather than RNAseq

## Required packages

- [optparse](https://github.com/trevorld/r-optparse)
- [tidyverse](https://www.tidyverse.org/)
- [biovisr](https://github.com/richysix/biovisr)
- [rnaseqtools](https://github.com/richysix/rnaseqtools)
