### Analysis of Enterovirus D68 (EV-D68) seroprevalence in England

This repo contains the data and code to run the sero-catalytic models presented in the paper by

Margarita Pons-Salort, Ben Lambert, Everlyn Kamau, Richard Pebody, Heli Harvala, Peter Simmonds, Nicholas C Grassly (2023) *Changes in transmission of Enterovirus D68 (EV-D68) in England inferred from seroprevalence data* eLife 12:e76609, available [here](https://doi.org/10.7554/eLife.76609) 

There are two models that include maternal antibodies and no seroreversion. The first model (`src/stan/model_ent68_mab_const.stan`) assumes the force of infection is constant over time, and the second model (`src/stan/model_ent68_mab_rw.stan`) allows it to change following a random walk of order one. Model fitting can be run from the script `src/R/model_ent68_main.R`.

