An EEGLAB-toolbox for finding [microstates](http://www.scholarpedia.org/article/EEG_microstates) in EEG signals.

Requires [EEGLAB](http://sccn.ucsd.edu/eeglab/) and [FastICA](http://research.ics.aalto.fi/ica/fastica/)

How to install:

Put all the files in the plugins/ folder, whereever EEGLAB is installed, and reload EEGLAB. 
Load a dataset.
A menu called "Microstates" should appear under Tools.

Currently implemented algorithms:

FastICA

K-means

Agglomerative/Hierarchical clustering

N-microstates + temporal smoothing (see [Pascual-Marqui et al. 1995](http://www.ncbi.nlm.nih.gov/pubmed/7622149))

Multi-garrote, a Variational Bayes version of the above algorithm.