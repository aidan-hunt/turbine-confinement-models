# blockage-correction
This repository contains MATLAB scripts for applying various blockage corrections to turbine performance data.

## Blockage corrector classes
`BWClosedChannel`: Closed-channel linear momentum and blockage correction as described by Barnsley and Wellicome (1990). Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135).
`HoulsbyOpenChannel`: Open-channel linear momentum and blockage correction as described by [Houlsby et al. (2008)](https://ora.ox.ac.uk/objects/uuid:5576d575-7bac-44b6-ac79-f698edcda40e). Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135).
`Dehtyriov Corrector`: A closed-channel two-scale blockage correction for turbine arrays as described by [Dehtyriov et al. (2023)](https://submissions.ewtec.org/proc-ewtec/article/view/366).
`WerleCorrector`: A closed-channel blockage correction as described by [Werle (2010)](https://arc.aiaa.org/doi/10.2514/1.44602). Follows the implementation of [Ross and Polagye (2020)](https://doi.org/10.1016/j.renene.2020.01.135).

## Helper classes
`BCBase`: superclass that defines several basic functions used by all of the blockage corrector classes.
`BCTranslator`: utility class that provides functions for translating cross-flow turbine data produced at UW into a format that the blockage corrector classes expect. If you plan to use these blockage corrector objects extensively, I recommend you write a similar sort of syntax translator code.
