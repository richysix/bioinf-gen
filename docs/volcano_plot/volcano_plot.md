# Create a volcano plot from a transcriptomics experiment

[Go to script](../../volcano_plot.R)

Script to produce a volcano plot from a RNA-seq experiment.  
It expects columns named `adjp` and `log2fc` in the input file.

``` bash
Rscript volcano_plot.R \
test_data/volcano-test-data.tsv \
volcano-basic.pdf
```

![Test volcano plot. It shows log10(Adjusted pvalue) plotted against
log2(Fold Change).](volcano-basic.png "Test volcano plot")

There are also options to label points that are above threshold for
Adjusted pvalue and/or log<sub>2</sub>(Fold Change).  
The `--labels` option turns on labelling.  
If `--labels` is set and neither of `--log2fc_threshold` or
`--pval_threshold` are set, the script will use default values of 2 and
1e-5 respectively.  
The `--log2fc_threshold` is a threshold on absolute log<sub>2</sub>(Fold
Change), so down-regulated genes will also be labelled.  
It is also possible to set just one of the thresholds.

``` bash
Rscript volcano_plot.R \
--labels --log2fc_threshold 5 \
--pval_threshold 1e-45 \
test_data/volcano-test-data.tsv \
volcano-labelled.pdf
```

![Volcano plot with some genes labelled. It shows log10(Adjusted pvalue)
plotted against log2(Fold
Change).](volcano-labelled.png "Labelled volcano plot")

The `--log2fc_or_pval` option changes the labelling behaviour to label
points that are over either threshold.  
The default behaviour is to only label points over both thresholds.

``` bash
Rscript volcano_plot.R \
--labels --log2fc_or_pval \
--log2fc_threshold 5 \
--pval_threshold 1e-45 \
test_data/volcano-test-data.tsv \
volcano-either-threshold.pdf
```

![Volcano plot with more genes labelled. It shows log10(Adjusted pvalue)
plotted against log2(Fold
Change).](volcano-either-threshold.png "Test volcano plot")

## Required packages

-   [optparse](https://github.com/trevorld/r-optparse)
-   [tidyverse](https://www.tidyverse.org/)
-   [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
-   [miscr](https://github.com/richysix/miscr)

## Optional packages

-   [ggrastr](https://cran.r-project.org/web/packages/ggrastr/vignettes/Raster_geoms.html) -
    Use to rasterise the points layer if installed
-   [svglite](https://svglite.r-lib.org/) - Used for svg output if
    installed
