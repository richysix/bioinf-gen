# Create a barchart of enriched GO terms

## go_barchart.R

[Go to script](../../go_barchart.R)

Script to produce a barchart from a file of GO enrichments. By default
the script expects columns named GO.ID, Term, FE, Set and up_down, but
these can be changed by setting options. The GO.IDs are plotted on the y
axis and the horizontal bars represent the Fold Enrichment (FE). The
bars are coloured by Set and depending on up_down are plotted to the
left or right.

There is an example file in the test_data directory of this repository.

``` bash
Rscript go_barchart.R test_data/test_data_go.tsv
```

![Bar chart of GO terms against Fold Enrichment. The bars are coloured
by GO domain](go-barchart.png "Default GO bar chart")

The column to use can be changed with the `--x_variable` option. In this
example only the top 20 terms (by the x variable) are plotted.

``` bash
Rscript go_barchart.R \
--x_variable log10p --x_axis_title="-log10[pvalue]" \
--fill_variable Set --top_terms 20 \
--output_file go_barchart_top20.svg \
test_data/test_data_go.tsv
```

![Bar chart of GO terms against -log10pvalue. The bars are coloured by
experiment](go-barchart-top20.png "GO bar chart of top 20 terms by -log10p")

## Required packages

- [tidyverse](https://www.tidyverse.org/)
- [grid](https://www.tidyverse.org/)
- [biovisr](https://github.com/richysix/biovisr)
- [miscr](https://github.com/richysix/miscr)
