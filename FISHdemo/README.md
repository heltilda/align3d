These files reproduce the results from the paper "Measuring chromosome conformation by fluorescence microscopy" by B. C. Ross et al. (preprint at https://www.biorxiv.org/content/10.1101/798512v1).

To reproduce results, install align3d together with the source files in this directory, using installations instructions in align3d/README.

Included here:

* align3d source files (.cpp, .cicada)
* analyzeReconstructions.m -- plots reconstruction performance relative to controls
* BinProbs.m -- plots p-value vs. frequency histograms
* compareLabelings.m -- compares performance of labeling 1 vs labeling 2 for experiments that attempted 2 labelings
* guessAllConformations.m (calls guessConformation.m) -- guesses conformations from mapping probabilities
* plotLRhistogram.m -- produces distance function histograms
* plotContour.m -- plots an inferred conformation
* plotDots.m -- plots spots in 3D space or labeled loci along a contour
* plotProbs.m -- plots locus-to-spot mappings
* processRecoloredExpt.cicada -- align3d script that infers mapping probabilities from locus/spot files
* setupRecoloredExpt.cicada -- align3d script that produces locus/spot files along with control spot files
* sketchGaussModel.m -- sketches a Gaussian chain model (Figure 1C)
