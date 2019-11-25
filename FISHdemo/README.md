These files reproduce the results from the paper "Measuring chromosome conformation by fluorescence microscopy" by B. C. Ross et al.

To reproduce results, install align3d together with the source files in this directory, using installations instructions in align3d/README.

Included here:

* align3d source files (.cpp, .cicada)
* analyzeReconstructions.m -- plots reconstruction performance relative to controls
* compareLabelings.m -- compares performance of labeling 1 vs labeling 2 for experiments that attempted 2 labelings
* guessAllConformations.m (calls guessConformation.m) -- guesses conformations from mapping probabilities
* plotContour.m -- plots an inferred conformation
* processRecoloredExpt.cicada -- align3d script that infers mapping probabilities from locus/spot files
* setupRecoloredExpt.cicada -- align3d script that produces locus/spot files along with control spot files
* sketchGaussModel.m -- sketches a Gaussian chain model (Figure 1C)
