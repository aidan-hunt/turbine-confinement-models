# linear-momentum
This repository contains MATLAB scripts for applying various confined-flow linear momentum actuator disk theory (LMADT) models to turbine performance data to predict aspects of the flow field around the turbine.
These LMADT models can be used as the basis for analytical blockage corrections, which predict turbine performance in unconfined flow from measured or simulated performance in confined flow.
This repository was developed in support of the following manuscripts/journal articles:
- Hunt et. al., "Experimental validation of a linear momentum and bluff-body model for high-blockage cross-flow turbine arrays", 2024. *Submitted to Journal of Fluid Mechanics*. [Preprint available on arXiv](https://doi.org/10.48550/arXiv.2408.16705).

## LMADT model classes
Each LMADT model is implemented as a MATLAB class:
- `BWClosedChannel`: Closed-channel linear momentum and blockage correction as described by Barnsley and Wellicome (1990). Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135).
- `HoulsbyOpenChannel`: Open-channel linear momentum and blockage correction as described by [Houlsby et al. (2008)](https://ora.ox.ac.uk/objects/uuid:5576d575-7bac-44b6-ac79-f698edcda40e). Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135). Note that, as the Froude number approaches 0, the results of `HoulsbyOpenChannel` are equivalent to that of `BWClosedChannel`.
- `DehtyriovCorrector`: A closed-channel two-scale linear-momentum-based blockage correction for turbine arrays as described by [Dehtyriov et al. (2023)](https://submissions.ewtec.org/proc-ewtec/article/view/366).

## Helper classes
- `BCBase`: superclass that defines several basic functions used by all of the LMADT or blockage corrector classes.
- `BCTranslator`: utility class that provides functions for translating cross-flow turbine data produced at the University of Washington into a format that the blockage corrector classes expect. If you plan to use these blockage corrector objects extensively, I recommend you write a similar sort of syntax translator code for your data format.

## Miscellaneous
- `WerleCorrector`: A closed-channel blockage correction as described by [Werle (2010)](https://arc.aiaa.org/doi/10.2514/1.44602). Note that this correction is NOT based on LMADT, but is simply provided for reference. Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135).
