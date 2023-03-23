### Analysis of Enterovirus D68 (EV-D68) seroprevalence in England

This repo contains the data and code to run the sero-catalytic models presented in the pre-print *Changes in transmission of Enterovirus D68 (EV-D68) in England inferred from seroprevalence data*.

There are two models that include maternal antibodies and no seroreversion. The first model (`src/stan/model_ent68_mab_const.rstan`) assumes the force of infection is constant over time, and the second model (`src/stan/model_ent68_mab_rw.rstan`) allows it to change following a random walk of order one. Model fitting can be run from the script `src/R/model_ent68_main.R`.

