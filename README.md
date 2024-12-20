#  mammal-SSD-code

Code and (some) data for Jones & Sheard 2023, https://doi.org/10.1098/rspb.2023.1211 , originally posted here for the benefit of the reviewers. The rest of the data can be found here: https://rs.figshare.com/collections/Supplementary_material_from_The_macroevolutionary_dynamics_of_mammalian_sexual_size_dimorphism_/6904490. 

Anyway, this is code to run the Bayesian phylogenetic mixed models, spread out into separate folders for the three logistic regressions and the several varieties of categorial regressions, as well as the BayesTraits 'Multistate' models. You'll need some additional functions for the categorical regressions; these are in the 'files needed to run categorical regressions' folder. If you're just looking for the main model presented in the manuscript, it's under 'models on the full dataset - no sociality' (a script that will run versions with and without the seasonality variables). Note that each of these models take several hours to run and that it's very easy to accidentally overwrite the output of previous models.

We also include a random 100 of the 10,000 trees in the Upham et al. 2019 (https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494) phylogeny, of which the MCMCglmm codes will trim as part of the initial data cleaning. The trees necessary to run the BayesTraits models are in that folder.
