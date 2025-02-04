# West_Nile_model_Guadeloupe
This code reproduces the following paper:

**Reconstructing the silent circulation of West Nile Virus in a Caribbean island during 15 years using sentinel serological data**

Celia Hamouche, Jennifer Pradel, Nonito Pagès, Véronique Chevalier, Sylvie Lecollinet, Jonathan Bastard *, Benoit Durand *

* **load_dataset.R**: imports and prepares the serological data
* **figures_Guadeloupe.R**: generates Figures 1 and 2, and Table 2
* **mosq_guadeloupe.R**: imports and prepares the entomological data, and uses it to fit the model (defined in mod_longitu_mosquito_Guad_2.BUGS) using MCMC (Step 1)
* **mod_estim_diag.R**: fits the serological models (defined in the four .BUGS scripts) to the data using MCMC (Step 2) and performs chains diagnostic
* **run_model.R**: runs the serological models using estimated parameters and generates Figure 3 of the paper
