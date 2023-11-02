# Produce a bubble plot from a GO enrichment

[Go to script](../../go_bubble_plot.R)

Script to produce a bubble plot using the output from Ian’s topgo
script. It expects files called “BP.sig.tsv”, “CC.sig.tsv” and
“MF.sig.tsv” in the working directory.

There is some test toy GO data in the test_data directory of this
repository. For example, to run the script with defaults

``` bash
cd test_data
Rscript ../go_bubble_plot.R
```

This will produce a bubble plot (go_bubble.pdf) with the top 5 terms by
pvalue labelled.

![Bubble plot of GO terms against -log10(pvalue). The points are
coloured by GO domain and the top 5 are
labelled](go-bubble-default.png "Default GO bubble plot")

To set a p value cut off for labelling

``` bash
cd test_data
Rscript ../go_bubble_plot.R --label_p_cutoff 1e-6
```

![Bubble plot of GO terms against -log10(pvalue). The points are
coloured by GO domain and points with pvalues below 1e-6 are
labelled](go-bubble-pval.png "GO bubble plot, terms with pvalue less than 1e-6 labelled")

Or to label specific terms. The GO IDs are used to specify which terms
to label, but the actual term descriptions are used as the labels.

``` bash
cd test_data
Rscript ../go_bubble_plot.R \
--labels "GO:0000001,GO:0000002,GO:0000003,GO:0000004,GO:0000005"
```

![Bubble plot of GO terms against -log10(pvalue). The points are
coloured by GO domain and the first 5 terms are
labelled](go-bubble-labels.png "GO bubble plot, with terms 1 to 5 labelled")

--no_labels will remove labels altogether

    ../go_bubble_plot.R --no_labels

## Required packages

- [tidyverse](https://www.tidyverse.org/)
- [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
- [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
- [biovisr](https://github.com/richysix/biovisr)
- [miscr](https://github.com/richysix/miscr)
