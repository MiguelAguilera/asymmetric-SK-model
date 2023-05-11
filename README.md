# A unifying framework for mean field theories of asymmetric kinetic Ising systems

Code reproducing the models used in the article Aguilera, M, Igarashi, M & Shimazaki H (2023). [Nonequilibrium thermodynamics of the asymmetric Sherrington-Kirkpatrick model](https://arxiv.org/abs/2205.09886).
<!--_Nature Communications_ 12:1197; [https://doi.org/10.1038/s41467-021-20890](https://doi.org/10.1038/s41467-021-20890).-->

## Abstract

Most natural systems operate far from equilibrium, displaying time-asymmetric, irreversible dynamics characterized by a positive entropy production while exchanging energy and matter with the environment. Although stochastic thermodynamics underpins the irreversible dynamics of small systems, the nonequilibrium thermodynamics of larger, more complex systems remains unexplored. Here, we investigate the asymmetric Sherrington-Kirkpatrick model with synchronous and asynchronous updates as a prototypical example of large-scale nonequilibrium processes. Using a path integral method, we calculate a generating functional over trajectories, obtaining exact solutions of the order parameters, path entropy, and steady-state entropy production of infinitely large networks. Entropy production peaks at critical order-disorder phase transitions, but is significantly larger for quasi-deterministic disordered dynamics. Consequently, entropy production can increase under distinct scenarios, requiring multiple thermodynamic quantities to describe the system accurately. These results contribute to developing an exact analytical theory of the nonequilibrium thermodynamics of large-scale physical and biological systems and their phase transitions.

## Description of the code

The main code reproduces the figures in the manuscript, which are saved in the 'img/' folder.

The data necessary to generate these figures is stored in the 'data/' folder, together with scripts for generating it. Note that the scripts 'generate-data-parallel-Glauber.py' and 'generate-data-sequential-Glauber.py' must be launched several times, which is computationally costy and was obtained parallelizing the code in a computer cluster.
