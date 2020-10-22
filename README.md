# bioinf-gen

General Bioinformatics scripts

## Scripts

* aggregate_zfa.R - aggregate ZFA terms by another column
* cluster_samples_rnaseq.R
* convert_detct_to_rnaseq_for_gsea.pl
* create_geneset_file.R
* create_rnaseq_rds.R
* deseq2-multiple_groups.R
* gene_expr_heatmap.R
* gene_lists_from_groups_cluego.pl
* get_msigdb_geneset.R
* go_bubble_plot.R
* graph_counts_by_group_facet.R - jittered and facetted count plot
* graph_counts_line.R
* histogram.R
* matrix_heatmap_plot.R
* merge_deseq_counts.pl - merge deseq counts from mutliple files
* mutmap_create_tsv.pl
* reshape-long_to_wide.R - reshapes long data to wide data
* volcano_plot.R
* xlsx_conditional_formatting.R

### go_bubble_plot.R

There is some test toy GO data in the test_data directory of this repository.
For example, run the script with defaults
```
cd test_data
../
../go_bubble_plot.R
```
This will produce a bubble plot (go_bubble.pdf) with the top 5 terms by pvalue labelled.
To set a p value cut off for labelling
```
../go_bubble_plot.R --label_p_cutoff 1e-6
```
Or to label specific terms
```
../go_bubble_plot.R \
--labels="GO:0000001,GO:0000002,GO:0000003,GO:0000004,GO:0000005"
```
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
