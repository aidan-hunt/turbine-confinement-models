# turbine-confinement-models
This repository contains MATLAB scripts for applying various confinement models to turbine performance data to predict aspects of the flow field around the turbine.
These models can be used as the basis for analytical blockage corrections, which predict turbine performance in unconfined flow from measured or simulated performance in confined flow.
This repository was developed in support of the following manuscripts/journal articles:
- Hunt et. al., "Experimental validation of a linear momentum and bluff-body model for high-blockage cross-flow turbine arrays". [Preprint available on arXiv](https://doi.org/10.48550/arXiv.2408.16705).
- Hunt et. al., "Performance characteristics and bluff-body modeling of high-blockage cross-flow turbine arrays with varying rotor geometry". [Preprint available on arXiv](https://doi.org/10.48550/arXiv.2410.19165).

## Confinement model classes
Each confinement model is implemented as a MATLAB class:
- `BWClosedChannel`: Closed-channel linear momentum and blockage correction developed by Barnsley and Wellicome (1990). Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135).
- `HoulsbyOpenChannel`: Open-channel linear momentum developed by [Houlsby et al. (2008)](https://ora.ox.ac.uk/objects/uuid:5576d575-7bac-44b6-ac79-f698edcda40e), with the blockage correction formulation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135). Follows the implementation of Ross and Polagye (2020). Note that, as the Froude number approaches 0, the results of `HoulsbyOpenChannel` are equivalent to that of `BWClosedChannel`.
- `NWTwoScale`: A closed-channel two-scale linear-momentum-based blockage correction for turbine arrays as described by [Nishino and Willden (2012)](https://doi.org/10.1017/jfm.2012.349), with the associated blockage correction developed by [Dehtyriov et al. (2023)](https://submissions.ewtec.org/proc-ewtec/article/view/366). Model implementation follows that of Dehtyriov et al (2023). NOTE: in this implementation, for consistency across induction factors we use an axial induction factor of *alpha = ut/Uinf* (where ut is the velocity at the turbine and Uinf is the freestream velocity), rather than *alpha = 1 - ut/Uinf* in Nishino and Willden (2012) or Dehtyriov et al (2023).
- `SteirosPotentialFlow`: A closed-channel confinement model and blockage correction based on two-dimensional potential flow theory developed by [Steiros et al. (2022)](https://doi.org/10.1017/jfm.2022.735).

## Helper classes
- `BCBase`: superclass that defines several basic functions used by all of the LMADT or blockage corrector classes.
- `BCTranslator`: utility class that provides functions for translating cross-flow turbine data produced at the University of Washington into a format that the blockage corrector classes expect. If you plan to use these blockage corrector objects extensively, I recommend you write a similar sort of syntax translator code for your data format.

## Miscellaneous
- `WerleCorrector`: A closed-channel blockage correction as described by [Werle (2010)](https://arc.aiaa.org/doi/10.2514/1.44602). Note that this correction is NOT based on LMADT, but is simply provided for reference. Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135).
