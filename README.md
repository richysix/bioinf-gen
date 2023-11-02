# bioinf-gen

General Bioinformatics scripts

## Scripts

Click on the links to go to the documentation for that script.

* aggregate_zfa.R - aggregate ZFA terms by another column
* [bubble_plot.R](docs/bubble_plot/bubble_plot.md) - Generic bubble plot
* cluster_samples_rnaseq.R
* convert_detct_to_rnaseq_for_gsea.pl
* create_geneset_file.R
* create_rnaseq_rds.R
* deseq2-multiple_groups.R
* [gene_expr_heatmap.R](docs/gene_expr_heatmap/gene_expr_heatmap.md) - Gene expression heatmap
* gene_lists_from_groups_cluego.pl
* get_msigdb_geneset.R
* [go_barchart.R](docs/go_barchart/go_barchart.md) - Produce a bar chart of GO results
* [go_bubble_plot.R](docs/go_bubble_plot/go_bubble_plot.md) - Produce a bubble plot from a topgo analysis
* [graph_counts_by_group_facet.R](docs/graph_counts_by_group_facet/graph_counts_by_group_facet.md) - jittered and facetted count plot
* graph_counts_line.R
* [gsea_to_genes.py](docs/gsea_to_genes/gsea_to_genes.md) - Get the genes behind GSEA enrichments
* histogram.R
* matrix_heatmap_plot.R
* merge_deseq_counts.pl - merge deseq counts from mutliple files
* mutmap_create_tsv.pl
* reshape-long_to_wide.R - reshapes long data to wide data
* [run_cluego.R](docs/run_cluego/run_cluego.md) - Run a Cytoscape ClueGO analysis from gene list(s)
* [upset-sig-genes.R](docs/upset-sig-genes/upset-sig-genes.md) - Simple UpSet plot from DESeq2 sig files with genes for each intersection
* [volcano_plot.R](docs/volcano_plot/volcano_plot.md) - Volcano plot
* xlsx_conditional_formatting.R

The plot scripts (gene_expr_heatmap.R, go_barchart.R, go_bubble_plot.R etc.) support outputting plot files as
pdf, png, ps and svg
