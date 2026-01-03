# RFBS
Bayesian phylogenetic model that extends the Dispersal-Extinction-Cladogenesis model to consider how anagenetic and cladogenetic events cause established and enabled biome affinities (or, more generally, other discrete realized versus fundamental niche states) to shift over evolutionary timescale. 

Directory includes:

- Sources files for RFBS in /src
- Scripts for performing analyses and generating plots from the RFBS paper. I included all scripts used to perform the anlyses and generate the figures included within the paper. However many of these are admittedly not very user-friendly. 
- Data file for Viburnum Analyses, including the phylogenies and data on observed estbalished affinities, included enabled affinities and excludes enabled affinities in /data/emp/viburnum. 


To make RFBS more user frinedly, I have also included a set of wrapper functions and an example script to "plug and play" given the users file name and phylogeny. An additional script is provided for users to generate ancestral affinitiy reconstruction figures as we did in the paper. MCMC chains can also be viewed using the program Tracer (downloaded here https://github.com/beast-dev/tracer/releases/tag/v1.7.2). 
