# bioinf-gen

General Bioinformatics scripts

## Scripts

* aggregate_zfa.R - aggregate ZFA terms by another column
* [bubble_plot.R](https://github.com/richysix/bioinf-gen#bubble_plotr) - Generic bubble plot
* cluster_samples_rnaseq.R
* convert_detct_to_rnaseq_for_gsea.pl
* create_geneset_file.R
* create_rnaseq_rds.R
* deseq2-multiple_groups.R
* [gene_expr_heatmap.R](https://github.com/richysix/bioinf-gen#gene_expr_heatmapr) - Gene expression heatmap
* gene_lists_from_groups_cluego.pl
* get_msigdb_geneset.R
* [go_barchart.R](https://github.com/richysix/bioinf-gen#go_barchartr) - Produce a bar chart of GO results
* [go_bubble_plot.R](https://github.com/richysix/bioinf-gen#go_bubble_plotr) - Produce a bubble plot from a topgo analysis
* [graph_counts_by_group_facet.R](https://github.com/richysix/bioinf-gen#graph_counts_by_group_facetr) - jittered and facetted count plot
* graph_counts_line.R
* gsea_to_genes.py - Get the genes behind GSEA enrichments
* histogram.R
* matrix_heatmap_plot.R
* merge_deseq_counts.pl - merge deseq counts from mutliple files
* mutmap_create_tsv.pl
* reshape-long_to_wide.R - reshapes long data to wide data
* [run_cluego.R](https://github.com/richysix/bioinf-gen#run_cluegor) - Run a Cytoscape ClueGO analysis from gene list(s)
* volcano_plot.R
* xlsx_conditional_formatting.R

The plot scripts (gene_expr_heatmap.R, go_barchart.R, go_bubble_plot.R) support outputting plot files as
pdf, png, ps and svg

### bubble_plot.R

Script to produce a bubble plot from continuous or categorical data. It expects 1 input file.
There is some test data in the test_data directory of this repository.
The script defaults to using columns named `x`, `y`, `fill` and `size` in the
data file. For example ...
```
cd test_data
../bubble_plot.R \
--output_file test_bubble_cont.pdf \
bubble_continuous.tsv
```

![Test bubble plot. It shows bubbles of different sizes and colours plotted
at random positions on the x and y axes.](test_data/test_bubble_cont.png "Test bubble plot")

Or for categorical data. This example also shows some of the possible options.
```
../bubble_plot.R --x_var Expt --y_var GO.ID \
--fill_var log10p --size_var Enrichment \
--y_labels Name --reverse_y \
--width 4 --height 4 \
--output_file test_bubble_cat.pdf \
go-3expts.tsv
```

![Test categorical bubble plot. It shows bubbles of different sizes and colours.
The x axis represents different experiments and the y axis represents enriched
Gene Ontology terms.](test_data/test_bubble_cat.png "Test Categorical bubble plot")

### gene_expr_heatmap.R

Script to produce a heatmap from RNAseq data. It expects a sample file and a count file (e.g. sig.tsv)
and the name of the output image file.

There is an example samples file and sig file in the test_data directory of this repository.
For example
```
../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
test_samples.tsv test_rnaseq_data.tsv test_heatmap.pdf
```

The available fill_palettes are those from [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html) and
[ColorBrewer](https://colorbrewer2.org/) via [scale_fill_distiller](https://ggplot2.tidyverse.org/reference/scale_brewer.html)

![Gene expression heatmap. Genes are displayed in rows with the samples in the
columns. Each box is coloured according to the expression of the gene/sample
combination](test_data/rnaseq_heatmap.png "RNAseq heatmap")

It is also possible to supply a list of gene ids to subset the heatmap to.

```
# create a test ids file
echo -e "ZFG005\nZFG006\nZFG009" > test_genes.txt

../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
--genes_file test_genes.txt \
--width 7 --height 2.5 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_subset.pdf
```

![Gene expression heatmap. Genes are displayed in rows with the samples in the
columns. Each box is coloured according to the expression of the gene/sample
combination. Only three of the original rows are shown](test_data/rnaseq_heatmap_subset.png "RNAseq heatmap subset")

Or only label the supplied genes
```
../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--sample_names --genes_to_label test_genes.txt \
--width 7 --height 2.5 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_subset.pdf
```

A file of sample metadata can also be supplied and will be plotted as a heatmap
under the expression heatmap.
```
../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
--sample_metadata_file test_samples_long.tsv \
--sample_metadata_ycol category \
--sample_metadata_fill value \
--relative_plot_sizes 9,2 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_with_sample_metadata.pdf
```

![Gene expression heatmap with sample metadata heatmap. Genes are displayed in
rows with the samples in the columns. Each box is coloured according to the
expression of the gene/sample combination. A second heatmap shows the metadata
associated with each sample](test_data/rnaseq_heatmap_with_sample_metadata.png "RNAseq heatmap with metadata")

The `--gene_tree` and `--sample_tree` options plot trees for the gene/sample clustering.

```
../gene_expr_heatmap.R --transform rlog \
--centre_and_scale --cluster both \
--colour_palette magma --cell_colour grey80 \
--gene_tree --sample_tree --gene_names --sample_names \
--sample_metadata_file test_samples_long.tsv \
--sample_metadata_ycol category --sample_metadata_fill value \
--relative_heights 1,3,1 --relative_widths 1,3 \
--width 480 --height 640 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_with_trees.png
```

![Gene expression heatmap with sample metadata heatmap and clustering trees. 
Genes are displayed in rows with the samples in the columns. Each box is 
coloured according to the expression of the gene/sample combination. 
A second heatmap shows the metadata associated with each sample and the 
clustering trees for the genes and samples are shown above and to the 
left of the main expression heatmap](test_data/test_heatmap_with_trees.png "RNAseq heatmap with sample metadata and clustering trees")

A file of gene metadata can also be supplied and will be plotted as a heatmap
to the right of the expression heatmap.

```
../gene_expr_heatmap.R --transform rlog \
--centre_and_scale --cluster both \
--colour_palette magma --cell_colour grey80 \
--gene_tree --sample_tree --gene_names --sample_names \
--sample_metadata_file test_samples_long.tsv \
--sample_metadata_ycol category --sample_metadata_fill value \
--gene_metadata_file test_gene_metadata.tsv \
--gene_metadata_xcol category --gene_metadata_fill value \
--relative_heights 1,3,1 --relative_widths 1,3,1 \
--width 640 --height 480 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_with_trees-sample_gene_metadata.png
```

![Gene expression heatmap with gene and sample metadata heatmaps and clustering trees. 
Genes are displayed in rows with the samples in the columns. Each box is coloured 
according to the expression of the gene/sample combination. Two more heatmaps show the 
metadata associated with each gene and sample and the clustering trees for the genes and 
samples are shown above and to the left of the main expression heatmap](test_data/test_heatmap_with_trees-sample_gene_metadata.png "RNAseq heatmap with gene and sample metadata and clustering trees")

### go_barchart.R

Script to produce a barchart from a file of GO enrichments. By default the
script expects columns named GO.ID, Term, FE, Set and up_down, but these can be
changed by setting options. The GO.IDs are plotted on the y axis and the
horizontal bars represent the Fold Enrichment (FE). The bars are coloured by
Set and depending os up_down are plotted to the left or right.

There is an example file in the test_data directory of this repository.

```
cd test_data
../go_barchart.R test_data_go.tsv
```

![Bar chart of GO terms against Fold Enrichment. The bars are coloured by GO
domain](test_data/go_barchart.png "Default GO bar chart")

The column to use can be changed with the `--x_variable` option. In this
example only the top 20 terms (by the x variable) are plotted.
```
../go_barchart.R --x_variable log10p --x_axis_title="-log10[pvalue]" \
--fill_variable Set --top_terms 20 \
--output_file go_barchart_top20.svg \
test_data_go.tsv
```

![Bar chart of GO terms against -log10pvalue. The bars are coloured by
experiment](test_data/go_barchart_top20.png "GO bar chart of top 20 terms by -log10p")

### go_bubble_plot.R

Script to produce a bubble plot using the output from Ian's topgo script.
It expects files called "BP.sig.tsv", "CC.sig.tsv" and "MF.sig.tsv" in the
working directory.

There is some test toy GO data in the test_data directory of this repository.
For example, run the script with defaults
```
cd test_data
../go_bubble_plot.R
```
This will produce a bubble plot (go_bubble.pdf) with the top 5 terms by pvalue labelled.

![Bubble plot of GO terms against -log10(pvalue). The points are coloured by GO
domain and the top 5 are labelled](test_data/go_bubble_plot_default.png "Default GO bubble plot")

To set a p value cut off for labelling
```
../go_bubble_plot.R --label_p_cutoff 1e-6
```

![Bubble plot of GO terms against -log10(pvalue). The points are coloured by GO
domain and points with pvalues below 1e-6 are labelled](test_data/go_bubble_plot_pval_threshold.png "GO bubble plot, terms with pvalue less than 1e-6 labelled")

Or to label specific terms. The GO IDs are used to specify which terms to label,
but the actual term descriptions are used as the labels.
```
../go_bubble_plot.R \
--labels="GO:0000001,GO:0000002,GO:0000003,GO:0000004,GO:0000005"
```

![Bubble plot of GO terms against -log10(pvalue). The points are coloured by GO
domain and the first 5 terms are labelled](test_data/go_bubble_plot_specific_labels.png "GO bubble plot, with terms 1 to 5 labelled")

--no_labels will remove labels altogether
```
../go_bubble_plot.R --no_labels
```

**Required packages**
* [tidyverse](https://www.tidyverse.org/)
* [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
* [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
* [biovisr](https://github.com/richysix/biovisr)
* [miscr](https://github.com/richysix/miscr)

### graph_counts_by_group_facet.R

This is a script to produce counts plots from RNA-seq (or DETCT) data.
The required arguments are a samples file and a sig file.
It expects the sample file to have a header with the column names.
It expects one of the columns to be called "sample".
e.g.
```
sample  condition
sample_1    wt
sample_2    mut
```

There is an example samples file and sig file in the test_data directory of
this repository.
The simplest way to run the script would be this:
```
cd test_data
../graph_counts_by_group_facet.R \
test_samples.tsv test_rnaseq_data.tsv
```

![Basic count plot showing the normalised counts for the wt condition as blue
circles and the mutant condition as orange circles](test_data/count_plot_basic.png "Basic count plot")

By default it tries to use a column named condition in the samples file as the
x variable. The default is to colour the points by condition as well.

The --x_variable option allows you to name a column to use as the x axis and
you can also specify columns to use for colours (--colour_variable) and
shapes (--shape_variable).
The --crossbar option plots a bar to show the mean/median of each group.
It is also possible to supply a file of Ensembl gene ids and the script will
only make plots for those gene ids in the data.

The example below uses almost all the available options
```
# create a test ids file
echo -e "ZFG005\nZFG006\nZFG009" > test_genes.txt

../graph_counts_by_group_facet.R \
--output_file test-condition-sex-treatment.pdf \
--genes_file test_genes.txt \
--x_variable condition \
--colour_variable condition \
--colour_palette wt=#0000ff,mut=#ff0000 \
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

![Count plot showing the normalised counts for the wt condition in blue and the
mutant condition in red. The points are split by the control or treated. Sex is
displayed as circle for female and squares for male](test_data/count_plot_everything.png "Count plot by condition by treatment")

The other options are:
* --no_jitter - removes the jitter from the points
* --detct - input data is DeTCT rather than RNAseq

**Required packages**
* [tidyverse](https://www.tidyverse.org/)
* [biovisr](https://github.com/richysix/biovisr)
* [rnaseqtools](https://github.com/richysix/rnaseqtools)

### gsea_to_genes.py

This is a script to take GSEA output files and return the genes that are
responsible for each enrichment.
There are test file in the test_data directory of this repository.
```
cd test_data
../gsea_to_genes.py --comparison test test_gsea_report.xls
```
Each line has the supplied comparison so that multiple of these output files
can be concatenated.
The `--genes_file` option allows the output to be limited to only genes in the
supplied list (e.g. sig genes).
```
../gsea_to_genes.py --genes_file gsea-genes.txt --comparison test test_gsea_report.xls
```
The default is to output to STDOUT, but an output filename can be given after
the input file name

### run_cluego.R

This script runs a standard ClueGO analysis from the supplied gene lists. It assumes that the gene list
has no header and is Ensembl gene ids in the first column. If more than one gene list is supplied,
the script will produce two images, one coloured by group (Enriched Term) and one coloured by cluster (Gene List origin).
At the moment the script runs the analysis and saves an image(s) and the output files, but I
can't find a way to save the analysis as a ClueGO session. This must be done manually in Cytoscape.

Cytoscape (>v3.6+) must be open and Cytoscape Apps 'yFiles Layout Algorithms' and 'ClueGO' must be installed.

*Example*
```
# Run with 1 gene list
run_cluego.R --verbose --analysis_name groups \
--output_image_file=cluego-groups.svg \
--output_basename=cluego-groups \
set1.sig.genes

# Run with 2 gene lists to see the overlap
run_cluego.R --verbose --analysis_name overlaps \
--output_image_file=overlap-cluego.svg \
--output_basename=overlap-cluego \
set1.sig.genes set2.sig.genes
```

I've had problems trying to run multiple analyses in sequence. If that happens you can try adding the `--destroy_network` option. 
That will destroy the network at the end of the script so you won't be able to save it as a ClueGO/Cytoscape network, but it means that you can produce a bunch of network pictures in one go.

**Required packages**
* [xml2](https://github.com/r-lib/xml2)
* [RJSONIO](https://cran.r-project.org/web/packages/RJSONIO/index.html)
* [httr](https://github.com/r-lib/httr)
* [biovisr](https://github.com/richysix/biovisr)
