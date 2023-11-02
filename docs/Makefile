md_files = ../test_data graph_counts_by_group_facet/graph_counts_by_group_facet.md \
gene_expr_heatmap/gene_expr_heatmap.md bubble_plot/bubble_plot.md \
go_bubble_plot/go_bubble_plot.md go_barchart/go_barchart.md \
volcano_plot/volcano_plot.md gsea_to_genes/gsea_to_genes.md \
run_cluego/run_cluego.md upset-sig-genes/upset-sig-genes.md \
kappa-scores-topgo/kappa-scores-topgo.md

# Documentation
docs: $(md_files)

../test_data: 
	mkdir ../test_data

graph_counts_files = graph_counts_by_group_facet/graph_counts_by_group_facet.md \
../test_data/test_samples.tsv ../test_data/test_rnaseq_data.tsv

$(graph_counts_files) : \
graph_counts_by_group_facet/graph_counts_by_group_facet.Rmd
	Rscript -e 'library("rmarkdown"); render("graph_counts_by_group_facet/graph_counts_by_group_facet.Rmd")'

gene_expr_heatmap/gene_expr_heatmap.md : gene_expr_heatmap/gene_expr_heatmap.Rmd \
../test_data/test_samples.tsv ../test_data/test_rnaseq_data.tsv
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

run_cluego/run_cluego.md: run_cluego/run_cluego.Rmd
	Rscript -e 'library("rmarkdown"); render("run_cluego/run_cluego.Rmd")'

upset-sig-genes/upset-sig-genes.md : upset-sig-genes/upset-sig-genes.Rmd
	Rscript -e 'library("rmarkdown"); render("upset-sig-genes/upset-sig-genes.Rmd")'

kappa-scores-topgo/kappa-scores-topgo.md : kappa-scores-topgo/kappa-scores-topgo.Rmd
	Rscript -e 'library("rmarkdown"); render("kappa-scores-topgo/kappa-scores-topgo.Rmd")'

.PHONY : cleanall cleandata cleanmd
cleanall : cleandata cleanmd

cleandata :
	-rm -r ../test_data

cleanmd : 
	find ./ -type f -name "*.md" | xargs rm
