# align3d
chromosome conformation reconstructions based on imaged loci

TO INSTALL:
* install the GNU Scientific Library
* install the GNU Multiple Precision Arithmetic Library
* download Cicada script (cicadascript.sourceforge.net)
* download the align3d files into the Cicada directory
* * overwrite Cicada's 'userfn.cpp', 'user.cicada' and 'Makefile' with align3d's versions of these files
* from the command line, navigate to the Cicada directory and type 'make align3d'
* to run from the command line, navigate to the Cicada directory and type './align3d'

TO TEST:
* from the align3d command prompt, type:  run("testA3D")

References:
* Ross BC, Wiggins PA. Measuring chromosome conformation with degenerate labels. Physical Review E. 2012 Jul 19;86(1):011918.
* Ross BC, Costello JC. Improved inference of chromosome conformation from images of labeled loci. F1000Research. 2018;7.
* Ross BC, Anaclerio F, Lorusso N, Ventura M, Costello JC. Measuring chromosome conformation by fluorescence microscopy (preprint). bioRxiv. https://www.biorxiv.org/content/10.1101/798512v1
