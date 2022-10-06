## bubble_plot.R

[Go to script](../../bubble_plot.R)

Script to produce a bubble plot from continuous or categorical data. It
expects 1 input file. There is some test data in the test_data directory
of this repository. The script defaults to using columns named `x`, `y`,
`fill` and `size` in the data file. For example …

``` bash
Rscript bubble_plot.R \
--output_file test_bubble_cont.pdf \
test_data/bubble_continuous.tsv
```

![Test bubble plot. It shows bubbles of different sizes and colours
plotted at random positions on the x and y
axes.](test_bubble_cont.png "Test bubble plot")

Or for categorical data. This example also shows some of the possible
options.

``` bash
Rscript bubble_plot.R --x_var Expt --y_var GO.ID \
--fill_var log10p --size_var Enrichment \
--y_labels Name --reverse_y \
--width 4 --height 4 \
--output_file test_bubble_cat.pdf \
test_data/go-3expts.tsv
```

![Test categorical bubble plot. It shows bubbles of different sizes and
colours. The x axis represents different experiments and the y axis
represents enriched Gene Ontology
terms.](test_bubble_cat.png "Test Categorical bubble plot")
