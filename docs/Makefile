md_files = graph_counts_by_group_facet/graph_counts_by_group_facet.md \
            gene_expr_heatmap/gene_expr_heatmap.md \
            bubble_plot/bubble_plot.md go_bubble_plot/go_bubble_plot.md \
            go_barchart/go_barchart.md \
            volcano_plot/volcano_plot.md \
            gsea_to_genes/gsea_to_genes.md

# Documentation
docs: $(md_files)

graph_counts_by_group_facet/graph_counts_by_group_facet.md: \
graph_counts_by_group_facet/graph_counts_by_group_facet.Rmd
	Rscript -e 'library("rmarkdown"); render("graph_counts_by_group_facet/graph_counts_by_group_facet.Rmd")'

gene_expr_heatmap/gene_expr_heatmap.md: gene_expr_heatmap/gene_expr_heatmap.Rmd
	Rscript -e 'library("rmarkdown"); render("gene_expr_heatmap/gene_expr_heatmap.Rmd")'

volcano_plot/volcano_plot.md: volcano_plot/volcano_plot.Rmd
	Rscript -e 'library("rmarkdown"); render("volcano_plot/volcano_plot.Rmd")'

bubble_plot/bubble_plot.md: bubble_plot/bubble_plot.Rmd
	Rscript -e 'library("rmarkdown"); render("bubble_plot/bubble_plot.Rmd")'

go_bubble_plot/go_bubble_plot.md: go_bubble_plot/go_bubble_plot.Rmd
	Rscript -e 'library("rmarkdown"); render("go_bubble_plot/go_bubble_plot.Rmd")'

go_barchart/go_barchart.md: go_barchart/go_barchart.Rmd
	Rscript -e 'library("rmarkdown"); render("go_barchart/go_barchart.Rmd")'

gsea_to_genes/gsea_to_genes.md: gsea_to_genes/gsea_to_genes.Rmd
	Rscript -e 'library("rmarkdown"); render("gsea_to_genes/gsea_to_genes.Rmd")'
