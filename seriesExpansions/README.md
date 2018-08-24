These files reproduce the results from the paper "Improved spatial alignmentfor inferring chromosome conformation" by B. C. Ross and J. Costello.

To reproduce results, install align3d using source files in THIS directory, using installations instructions in align3d/README.
Conformations were generated using:  github.com/heltilda/wormulator

Included here:
* align3d source files (.cpp, .cicada)
* exploreSeriesParams.cicada -- run this from wormulator to generate simulated conformations, or call from align3d to process them
* makeWormChains.cicada -- called by exploreSeriesParams.cicada when generating conformations
* plotDots.m -- plots simulated conformations with labelings (e.g. plots in Fig 3A)
* plotProbs.m -- plots locus-to-spot mapping probabilities (Fig 3B)
* plotSeriesComparison.m -- produces Figs 5, 6, S2, S3 of paper
* makePercentilePlot.m -- called by plotSeriesComparison.m
* makeCostFunctionHistograms.m -- produces Fig 4 of paper
* plotConformationError.m -- produces Fig 7 of paper
* BinProbs.m -- produces Fig S4 of paper
