# CFI-MI MATLAB toolbox

This implements the simulation of a new proposed measure by G. Mijatovic, T. Loncar Turukalo, N. Bozanic and L. Faes: 
"A measure of concurrent neural firing activity based on mutual information", 2020.

The CFI_MI index estimates the degree of concomitant firing between two neural units based on a modified form of mutual information (MI) 
applied on a coarse, binary representations of firing activity (states of neural quiescence "0" and states of firing "1" unfolded in time).


The CFI-MI toolbox contains next functions:

1. binary_representation.m: binary-state representation of the spiking activity (see REF0);
2. function_CFI_MI.m: computing the concurrent firing index based on modified mutual information, CFI_MI index;
3. spiSeMe_surrogate_jodi.m: assessment of statistical significance of the CFI_MI index index based on surrogate data analysis (this function is part of the SpiSeMe package, see REF2);

and script:

4. demo.m: to demonstrate the estimaton of the CFI_MI index between two selected cells from a network of 1000 randomly coupled spiking neurons (see REF1); supported with assessment of its statistical significance (see REF2) and the corresponding visualization of spike trains and binary profiles.

REF0: Mijatovic G, Loncar-Turukalo T, Procyk E, Bajic D (2018): "A novel approachto probabilistic characterisation of neural firing patterns", Journal of neuroscience methods 305:67–81

REF1: Izhikevich EM (2003): "Simple model of spiking neurons", IEEE Transactions onneural networks 14(6):1569–1572

REF2: Ricci L, Castelluzzo M, Minati L, Perinelli A  (2019): "Generation  of  surro-gate event sequences via joint distribution of successive inter-event intervals", Chaos: An Interdisciplinary Journal of Nonlinear Science 29(12):121102

