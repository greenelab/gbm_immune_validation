# Investigating immune cell profiles across Glioblastoma subtypes

**Gregory Way and Casey Greene**
**University of Pennsylvania**

## Summary

The amount and type of immune cell infiltration into tumors is an important
determinant of disease progression and survival. Cancer subtypes can often be
distinguished by the rate of immune cell infiltration since this tends to be one
of the more dominant observable signatures in gene expression data. Current
methods for directly observing immune cell profiles in a given population
include laboriously quantifying immune cell proportions by flow cytometry or
other technically challenging cell labeling techniques. Therefore, deconvolution
methods are being developed to automatically extract immune cell proportions
from full tumor gene expression data.

We used ssGSEA ([Barbie et al. 2009](http:/doi.org/10.1038/nature08460 "ssGSEA"))
to deconvolute immune cell signatures from glioblastoma multiforme (GBM) tumors
from The Cancer Genome Atlas. Briefly, ssGSEA is a simple rank based test that
evaluates the empirical cumulative distribution function of input gene sets
compared to the eCDF of the remaining genes.

We used LM22.txt as defined by
[Newman et al. 2015](http://doi.org/10.1038/nmeth.3337 "CIBERSORT") as input
genelists to ssGSEA.

## Reproducibility

Our end to end analysis from downloading data to generating publication ready
figures is provided in this github repository. We implement an automatic 
[reproducible workflow using continuous analysis](http://doi.org/10.1101/056473 "Beaulieu-Jones and Greene 2016")
to ensure a stable compute environment and consistent reproducibility.

We use the ssGSEA implementation available on
[bioconductor](https://bioconductor.org/packages/release/bioc/html/GSVA.html)
(Guinney and Castelo 2016).

```bash
# To reproduce the pipeline independently simply run:
bash run_pipeline.sh
```

For exact instructions on how to reproduce our analysis see `run_pipeline.sh`.

## Contact

For all code related questions, bug reporting, or feature requests please file a
[GitHub issue](https://github.com/greenelab/GBM_immune_profiles/issues "file an issue")

## Dependencies

All code was performed in R version 3.2.3 and packages were managed by the
[checkpoint](https://cran.r-project.org/web/packages/checkpoint/index.html "Microsoft checkpoint r package")
package set at "2016-08-16" (see `install.R` for more details)

