# Core code for QTL inference =======================================
This folder contains the scripts used for QTL inference from genotype and 
phenotype data.

step_QTL.py takes an optional starting QTL model, phenotype and genotype
data, and a number of QTLs to add to the model. It does it in a stepwise
fashion where at each step it finds the next best new QTL to add to the model,
and reoptimizes all QTL positions'.

step_epistasis.py does the same thing, but for pairwise epistatic terms.
It requires an additive model to start.
