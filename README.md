# RFBS
Bayesian phylogenetic model that extends the Dispersal-Extinction-Cladogenesis model to consider how anagenetic and cladogenetic events cause established and enabled biome affinities (or, more generally, other discrete realized versus fundamental niche states) to shift over evolutionary timescale. 

Directory includes:

- Sources files for RFBS in /src
- Scripts for performing analyses and generating plots from the RFBS paper. I included all scripts used to perform the anlyses and generate the figures included within the paper.  
- Data file for Viburnum Analyses, including the phylogenies and data on observed estbalished affinities, included enabled affinities and excludes enabled affinities in /data/emp/viburnum. 


Mny of these are admittedly not very user-friendly, RFBS originally accepted multiple data files for "standard" observed affintiies, excluded and included enabled affinities. To make RFBS more user frinedly and generally-applicable, I have included a set of wrapper functions and an example script to "plug and play". The user only needs to provide the phylogeny and a single affinity data file – a csv where the first column is species names that match the phylogenetic tip labels, and the subsequent columns score biome affinities (2=Established Affinities, 1=an explicitly included enabled affinity, 0=and explicitly excluded enabled affinity, and NA=an unobserved affinity that could either be an enabled affinity or a non-affinity). An additional script is provided for users to generate ancestral affinitiy reconstruction figures as we did in the paper. MCMC chains can also be viewed using the program Tracer (downloaded here https://github.com/beast-dev/tracer/releases/tag/v1.7.2). 
