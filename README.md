# motoneuron-matlab
This repository encompasses the source code for three motoneuron models with active dendrites.

# Description
Experimental results have pointed out that the dendritic tree of spinal motoneurons is not a passive resistive-capacitive network. Instead, voltage-gated ionic conductances are present along the dendritic tree accounting for many nonlinear behaviors, such as bistability (or a selfsustained firing following a brief excitatory input) and plateau potentials. Our aim is to present three new motoneuron models with a dendritic L-type calcium conductance so that a persistent inward current (PIC) could be generated, contributing to the genesis of phenomena observed in decerebrate cat motoneurons. The models are of types S, FR, and FF with two compartments each, one for the soma and another for the dendritic tree. The geometric and electrotonic parameters were based on <a href="http://dx.doi.org/10.1007/s10827-008-0092-8" target="_blank">a previous development</a> from our laboratory in which passive-dendrite motoneuron models were developed. Similarly to the passive-dendrite models, the soma compartment comprises a sodium channel and a fast potassium channel, both responsible for the genesis of action potential, along with a slow potassium channel yielding the afterhyperpolarization (AHP). In our new models, an L-type calcium channel was implemented in the denditic compartment and its parameters were chosen according the properties of PIC and the firing patterns of decerebrate cat motoneurons. The models were coded in <a href="http://www.mathworks.com" target="_blank">Matlab®</a>. We encourage you to use or to improve these models in accordance with the distribution license, but without forgetting to cite the source.

# Publications
1. Elias L.A. and Kohn A.F. (2013) Individual and collective properties of computationally efficient motoneuron models of types S and F with active dendrites. **Neurocomputing**, 99: 521-533. <a href="http://dx.doi.org/10.1016/j.neucom.2012.06.038" target="_blank">[doi]</a>
2. Elias L.A., Chaud V.M., and Kohn A.F. (2012) Models of passive and active dendrite motoneuron pools and their differences in muscle force control. **Journal of Computational Neuroscience**, 33(3): 515-531. <a href="http://dx.doi.org/10.1007/s10827-012-0398-4" target="_blank">[doi]</a>

# Funding
This project was fully funded by <a href="http://www.fapesp.br/en/" target="_blank">FAPESP</a> (project no. <a href="http://www.bv.fapesp.br/en/bolsas/109027/computer-simulation-of-motoneuron-behavior-with-passive-and-active-dendrites-in-response-to-synaptic/" target="_blank">2007/05898-4</a>)

# License
This software is licensed under the [GNU License v3](LICENSE).