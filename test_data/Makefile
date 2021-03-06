outputs: test_heatmap.pdf test_heatmap_subset.pdf \
test_heatmap_with_metadata.pdf test_heatmap_with_trees.pdf \
test_heatmap_with_trees-sample_gene_metadata.png \
go_bubble.pdf go_bubble_pval.pdf go_bubble_specific.pdf go_bubble_no_labels.pdf \
counts.pdf test-condition-sex-treatment.pdf \
go_barchart.pdf go_barchart_top20.svg \
test-gsea_to_genes.out.tsv test-gsea_to_genes.sig.out.tsv

# test data
BP.sig.tsv CC.sig.tsv MF.sig.tsv \
test_samples.tsv test_samples_long.tsv \
test_gene_metadata.tsv test_rnaseq_data.tsv \
test_data_go.tsv test_gsea_report.xls \
Term001.xls Term002.xls Term003.xls Term004.xls Term005.xls gsea-genes.txt \
go-3expts.tsv bubble_continuous.tsv: test_data.R
	./test_data.R

# bubble plots
test_bubble_cont.pdf: ../bubble_plot.R go-3expts.tsv
	../bubble_plot.R \
	--output_file test_bubble_cont.pdf \
	bubble_continuous.tsv

test_bubble_cat.pdf: ../bubble_plot.R go-3expts.tsv
	../bubble_plot.R --x_var Expt --y_var GO.ID \
	--fill_var log10p --size_var Enrichment \
	--y_labels Name --reverse_y \
	--width 4 --height 4 \
	--output_file test_bubble_cat.pdf \
	go-3expts.tsv


# Gene expression heatmaps
test_heatmap.pdf: ../gene_expr_heatmap.R test_samples.tsv test_rnaseq_data.tsv
	../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
test_samples.tsv test_rnaseq_data.tsv test_heatmap.pdf

test_genes.txt:
	echo -e "ZFG005\nZFG006\nZFG009" > test_genes.txt

test_heatmap_subset.pdf: ../gene_expr_heatmap.R test_genes.txt test_samples.tsv test_rnaseq_data.tsv
	../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
--genes_file test_genes.txt \
--width 7 --height 2.5 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_subset.pdf

test_heatmap_with_metadata.pdf: ../gene_expr_heatmap.R test_samples_long.tsv test_samples.tsv test_rnaseq_data.tsv
	../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster rows \
--colour_palette magma --cell_colour grey80 \
--gene_names --sample_names \
--sample_metadata_file test_samples_long.tsv \
--sample_metadata_ycol category \
--sample_metadata_fill value \
--relative_heights 0,9,2 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_with_metadata.pdf

test_heatmap_with_trees.pdf: ../gene_expr_heatmap.R test_samples_long.tsv test_samples.tsv test_rnaseq_data.tsv
	../gene_expr_heatmap.R --transform rlog \
--center_and_scale --cluster both \
--colour_palette magma --cell_colour grey80 \
--gene_tree --sample_tree \
--gene_names --sample_names \
--sample_metadata_file test_samples_long.tsv \
--sample_metadata_ycol category \
--sample_metadata_fill value \
--relative_heights 1,3,1 \
--relative_widths 1,3 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_with_trees.pdf

test_heatmap_with_trees-sample_gene_metadata.png: ../gene_expr_heatmap.R \
test_samples.tsv test_rnaseq_data.tsv
	../gene_expr_heatmap.R --transform rlog \
--centre_and_scale --cluster both \
--colour_palette magma --cell_colour grey80 \
--gene_tree --sample_tree --gene_names --sample_names \
--sample_metadata_file test_samples_long.tsv \
--sample_metadata_ycol category --sample_metadata_fill value \
--gene_metadata_file test_gene_metadata.tsv \
--gene_metadata_xcol value --gene_metadata_fill category \
--relative_heights 1,3,1 --relative_widths 1,3,1 \
--width 640 --height 480 \
test_samples.tsv test_rnaseq_data.tsv test_heatmap_with_trees-sample_gene_metadata.png

# GO bubble plots
go_bubble.pdf: ../go_bubble_plot.R BP.sig.tsv CC.sig.tsv MF.sig.tsv
	../go_bubble_plot.R

go_bubble_pval.pdf: ../go_bubble_plot.R BP.sig.tsv CC.sig.tsv MF.sig.tsv
	../go_bubble_plot.R --label_p_cutoff 1e-6 --output_file go_bubble_pval.pdf

go_bubble_specific.pdf: ../go_bubble_plot.R BP.sig.tsv CC.sig.tsv MF.sig.tsv
	../go_bubble_plot.R --output_file go_bubble_specific.pdf \
--labels="GO:0000001,GO:0000002,GO:0000003,GO:0000004,GO:0000005"

go_bubble_no_labels.pdf: ../go_bubble_plot.R BP.sig.tsv CC.sig.tsv MF.sig.tsv
	../go_bubble_plot.R --no_labels --output_file go_bubble_no_labels.pdf

# RNAseq count plots
counts.pdf: ../graph_counts_by_group_facet.R test_samples.tsv test_rnaseq_data.tsv
	../graph_counts_by_group_facet.R test_samples.tsv test_rnaseq_data.tsv

test-condition-sex-treatment.pdf: ../graph_counts_by_group_facet.R test_genes.txt test_samples.tsv test_rnaseq_data.tsv
	../graph_counts_by_group_facet.R \
--output_file test-condition-sex-treatment.pdf \
--genes_file test_genes.txt \
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

# GO barcharts
go_barchart.pdf: ../go_barchart.R test_data_go.tsv
	../go_barchart.R test_data_go.tsv
 
go_barchart_top20.svg: ../go_barchart.R test_data_go.tsv
	../go_barchart.R --x_variable log10p --x_axis_title="-log10[pvalue]" \
--fill_variable Set --top_terms 20 \
--output_file go_barchart_top20.svg \
test_data_go.tsv

test-gsea_to_genes.out.tsv: ../gsea_to_genes.py test_gsea_report.xls \
Term001.xls Term002.xls Term003.xls Term004.xls Term005.xls
	../gsea_to_genes.py --comparison test test_gsea_report.xls > test-gsea_to_genes.out.tsv

test-gsea_to_genes.sig.out.tsv: ../gsea_to_genes.py test_gsea_report.xls \
Term001.xls Term002.xls Term003.xls Term004.xls Term005.xls
	../gsea_to_genes.py --genes_file gsea-genes.txt --comparison test test_gsea_report.xls > test-gsea_to_genes.sig.out.tsv

.PHONY : clean
clean :
	-rm *.pdf *.svg

.PHONY : clean_all
clean_all :
	-rm *.tsv *.xls *.rda *.pdf *.svg
