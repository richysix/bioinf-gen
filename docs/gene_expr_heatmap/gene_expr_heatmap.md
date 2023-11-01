# Create gene expression heatmaps from a transcriptomics experiment

## gene_expr_heatmap.R

Script to produce a heatmap from RNAseq data. It expects a sample file
and a count file (e.g. sig.tsv) and the name of the output image file.

There is an example samples file and sig file in the test_data directory
of this repository. For a basic heatmap do:

``` bash
Rscript gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv test_data/test_heatmap.pdf
```

The available fill_palettes are those from
[viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
and [ColorBrewer](https://colorbrewer2.org/) via
[scale_fill_distiller](https://ggplot2.tidyverse.org/reference/scale_brewer.html)

![Gene expression heatmap. Genes are displayed in rows with the samples
in the columns. Each box is coloured according to the expression of the
gene/sample combination](rnaseq-heatmap.png "RNAseq heatmap")

It is also possible to supply a list of gene ids to subset the heatmap
to.

``` bash
# create a test ids file
echo -e "GeneID\tName
ENSTEST005\tgene-5
ENSTEST006\tgene-6
ENSTEST009\tgene-9" > test_data/test_genes.txt

Rscript gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
--genes_file test_data/test_genes.txt \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv test_data/test_heatmap_subset.pdf
```

![Gene expression heatmap. Genes are displayed in rows with the samples
in the columns. Each box is coloured according to the expression of the
gene/sample combination. Only three of the original rows are
shown](rnaseq-heatmap-subset.png "RNAseq heatmap subset")

Or only label the supplied genes

``` bash
Rscript gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--sample_names \
--genes_to_label test_data/test_genes.txt \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv test_data/test_heatmap_subset_labels.pdf
```

![Gene expression heatmap. Genes are displayed in rows with the samples
in the columns. Each box is coloured according to the expression of the
gene/sample combination. Only three of the rows are
labelled](rnaseq-heatmap-subset-labels.png "RNAseq heatmap with a subset of genes labelled")

A file of sample metadata can also be supplied and will be plotted as a
heatmap under the expression heatmap.

``` bash
Rscript gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
--sample_metadata_file test_data/test_samples_long.tsv \
--sample_metadata_ycol category \
--sample_metadata_fill value \
--relative_heights "0,8,2" \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv test_data/test_heatmap_with_sample_metadata.pdf
```

![Gene expression heatmap with sample metadata heatmap. Genes are
displayed in rows with the samples in the columns. Each box is coloured
according to the expression of the gene/sample combination. A second
heatmap shows the metadata associated with each
sample](rnaseq-heatmap-sample-metadata.png "RNAseq heatmap with metadata")

The `--gene_tree` and `--sample_tree` options plot trees for the
gene/sample clustering.

``` bash
Rscript gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster both \
--colour_palette magma --cell_colour grey80 \
--gene_tree --sample_tree --gene_names --sample_names \
--sample_metadata_file test_data/test_samples_long.tsv \
--sample_metadata_ycol category \
--sample_metadata_fill value \
--relative_heights 1,3,1 --relative_widths 1,3 \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv test_data/test_heatmap_with_trees.pdf
```

![Gene expression heatmap with sample metadata heatmap and clustering
trees. Genes are displayed in rows with the samples in the columns. Each
box is coloured according to the expression of the gene/sample
combination. A second heatmap shows the metadata associated with each
sample and the clustering trees for the genes and samples are shown
above and to the left of the main expression
heatmap](rnaseq-heatmap-trees.png "RNAseq heatmap with sample metadata and clustering trees")

A file of gene metadata can also be supplied and will be plotted as a
heatmap to the right of the expression heatmap.

``` bash
Rscript gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster both \
--colour_palette magma --cell_colour grey80 \
--gene_tree --sample_tree --gene_names --sample_names \
--sample_metadata_file test_data/test_samples_long.tsv \
--sample_metadata_ycol category \
--sample_metadata_fill value \
--gene_metadata_file test_data/test_gene_metadata.tsv \
--gene_metadata_xcol category --gene_metadata_fill value \
--gene_metadata_fill_palette cividis \
--relative_heights 1,3,1 --relative_widths 1,3,1 \
--width 640 --height 480 \
test_data/test_samples.tsv test_data/test_rnaseq_data.tsv test_heatmap_with_gene-metadata.pdf
```

![Gene expression heatmap with gene and sample metadata heatmaps and
clustering trees. Genes are displayed in rows with the samples in the
columns. Each box is coloured according to the expression of the
gene/sample combination. Two more heatmaps show the metadata associated
with each gene and sample and the clustering trees for the genes and
samples are shown above and to the left of the main expression
heatmap](rnaseq-heatmap-gene-metadata.png "RNAseq heatmap with gene and sample metadata and clustering trees")

## Required packages

- [tidyverse](https://www.tidyverse.org/)
- [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
- [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [patchwork](https://patchwork.data-imaginist.com/)
- [biovisr](https://github.com/richysix/biovisr)
- [rnaseqtools](https://github.com/richysix/rnaseqtools)
