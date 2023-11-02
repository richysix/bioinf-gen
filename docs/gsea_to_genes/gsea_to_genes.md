# Get genes for GSEA enrichments

## gsea_to_genes.py

[Go to script](../../gsea_to_genes.py)

This is a script to take GSEA output files and return the genes that are
responsible for each enrichment. There are test files in the test_data
directory of this repository.

``` bash
cd test_data
python ../gsea_to_genes.py --comparison test test_gsea_report.xls
```

Each line has the supplied comparison so that multiple of these output
files can be concatenated.

The `--genes_file` option allows the output to be limited to only genes
in the supplied list (e.g. sig genes).

``` bash
cd test_data
python ../gsea_to_genes.py --genes_file gsea-genes.txt --comparison test test_gsea_report.xls
```

The default is to output to STDOUT, but an output filename can be given
after the input file name

If running from a directory other than the one containing the GSEA
output files the `--base_dir` option should be set.

``` bash
python gsea_to_genes.py --genes_file test_data/gsea-genes.txt --comparison test \
--base_dir test_data test_data/test_gsea_report.xls
```
