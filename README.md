
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Compendium of R code and data for “Estimation of the Degree of Decomposition of Peat and Past Net Primary Production from Mid-Infrared Spectra”

This repository contains the data and code for our manuscript:

> Teickner, Henning, Julien Arsenault, Mariusz Gałka, and Klaus-Holger
> Knorr. 2025. “Estimation of the Degree of Decomposition of Peat and
> Past Net Primary Production from Mid-Infrared Spectra.” (unpublished).

### How to cite

Please cite this compendium as:

> Teickner, Henning, Julien Arsenault, Mariusz Gałka, and Klaus-Holger
> Knorr, (2025). Compendium of R code and data for “Estimation of the
> Degree of Decomposition of Peat and Past Net Primary Production from
> Mid-Infrared Spectra.” Accessed 26 Sep 2025.
> <https://github.com/henningte/eb1149>

### How to use

Instructions how to set up the Docker containers to reproduce the
computations are available from the Dockerfile. The Dockerfile also
provides instructions to run the
[`targets`](https://github.com/ropensci/targets) workflow to reproduce
all computations and eventually the manuscript and supporting
information.

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html)

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse. See the sources section for licenses for
data derived from external sources, their licenses, and how to give
credit to the original author(s) and the source.

### Sources

Data in folder `data/raw_data` are derived from different sources. To
use these data and give credit to data authors, please follow the
following information:

- [:file_folder: pmird](data/raw_data/pmird): This is a template folder
  in which the folder `pmird_prepared_data` from the pmird database
  needs to be stored to reproduce the computations. This folder is
  available from <https://doi.org/10.5281/zenodo.17092587>.

- [:file_folder:
  ecy2462-sup-0005-npp_moss.csv](data/raw_data/ecy2462-sup-0005-npp_moss.csv):
  This is one of the csv files contained in the Peatland Decomposition
  and Productivity Database (Bona et al. 2018). The reproduction is a
  copy of an official work that is published by Natural Resources Canada
  (NRCan) and has not been produced in affiliation with, or with the
  endorsement of, NRCan. Commercial reproduction and distribution are
  prohibited except with written permission from NRCan. For more
  information, contact NRCan at
  <copyright.droitdauteur@nrcan-rncan.gc.ca>.

  The full copyright statement is: © Her Majesty the Queen in Right of
  Canada, 2018. Information contained in this publication or product may
  be reproduced, in part or in whole, and by any means, for personal or
  public noncommercial purposes, without charge or further permission,
  unless otherwise specified. You are asked to exercise due diligence in
  ensuring the accuracy of the materials reproduced; indicate the
  complete title of the materials reproduced, and the name of the author
  organization; and indicate that the reproduction is a copy of an
  official work that is published by Natural Resources Canada (NRCan)
  and that the reproduction has not been produced in affiliation with,
  or with the endorsement of, NRCan. Commercial reproduction and
  distribution are prohibited except with written permission from NRCan.
  For more information, contact NRCan at
  <copyright.droitdauteur@nrcan-rncan.gc.ca>.

  See Bona et al. (2018)
  (<https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2462>)
  for details.

- [:file_folder:
  Bengtsson_etal_sph_holarctic_growth.csv](data/raw_data/Bengtsson_etal_sph_holarctic_growth.csv):
  csv file from Bengtsson et al. (2020) (these data are described in
  more detail in Bengtsson et al. (2021)). License:
  [CC-0](https://creativecommons.org/publicdomain/zero/1.0/). See
  <http://datadryad.org/stash/dataset/doi:10.5061/dryad.1ns1rn8rm> for
  details.

- [:file_folder: d88_2.rds, d89_2.rds, d90_1.rds](data/raw_data): These
  are average testate amoebae-reconstructed water table depth for the
  peat cores analyzed in our manuscript, digitized from figures in the
  original publications (Gałka, Diaconu, et al. 2022; Gałka, Hölzer, et
  al. 2022; Diaconu et al. 2020). License:
  [CC-0](https://creativecommons.org/publicdomain/zero/1.0/)

- [:file_folder:
  pangaea_metadata_reuter2019.rds](data/raw_data/pangaea_metadata_reuter2019.rds):
  rds file with the metadata from Reuter et al. (2019), extracted from
  the Pangaea repository using the `pangaear` R package (Chamberlain et
  al. 2025). License:
  [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/). For
  details, see <https://doi.pangaea.de/10.1594/PANGAEA.902181>.

- [:file_folder: eb1064](data/raw_data/eb1064): Folder with preprocessed
  data from Arsenault et al. (2024a) (the contents of the folder were
  not directly derived from Arsenault et al. (2024a), but processed from
  the same original raw data files) and containing in addition the
  mid-infrared spectra for the samples in Arsenault et al. (2024a) (for
  an overview on the data, see Arsenault et al. (2024b)).

### Contributions

We welcome contributions from everyone. Please note that the eb1149
project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

### Funding

This study was funded by the Deutsche Forschungsgemeinschaft (DFG,
German Research Foundation) grant no. KN 929/23-1 to Klaus-Holger Knorr
and grant no. PE 1632/18-1 to Edzer Pebesma.

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Arsenault.2024b" class="csl-entry">

Arsenault, Julien, Julie Talbot, Tim R. Moore, Klaus-Holger Knorr,
Henning Teickner, and Jean-François Lapierre. 2024a. “Patterns and
Drivers of Organic Matter Decomposition in Peatland Open-Water Pools.”
Zenodo. <https://doi.org/10.5281/ZENODO.10581235>.

</div>

<div id="ref-Arsenault.2024a" class="csl-entry">

———. 2024b. “Patterns and Drivers of Organic Matter Decomposition in
Peatland Open-Water Pools.” *Biogeosciences* 21 (15): 3491–3507.
<https://doi.org/10.5194/bg-21-3491-2024>.

</div>

<div id="ref-Bengtsson.2021" class="csl-entry">

Bengtsson, Fia, Håkan Rydin, Jennifer L. Baltzer, Luca Bragazza,
Zhao-Jun Bu, Simon J. M. Caporn, Ellen Dorrepaal, et al. 2021.
“Environmental Drivers of *Sphagnum* Growth in Peatlands Across the
Holarctic Region.” Edited by Rien Aerts. *Journal of Ecology* 109 (1):
417–31. <https://doi.org/10.1111/1365-2745.13499>.

</div>

<div id="ref-Bengtsson.2020" class="csl-entry">

Bengtsson, Fia, Håkan Rydin, Jennifer Baltzer, Luca Bragazza, Zhao-Jun
Bu, Simon Caporn, Ellen Dorrepaal, et al. 2020. “Environmental Drivers
of *Sphagnum* Growth in Peatlands Across the Holarctic Region.” Dryad.
<https://doi.org/10.5061/DRYAD.1NS1RN8RM>.

</div>

<div id="ref-Bona.2018" class="csl-entry">

Bona, Kelly Ann, Arlene Hilger, Magdalena Burgess, Nicole Wozney, and
Cindy Shaw. 2018. “A Peatland Productivity and Decomposition Parameter
Database.” *Ecology* 99 (10): 2406–6.
<https://doi.org/10.1002/ecy.2462>.

</div>

<div id="ref-Chamberlain.2025" class="csl-entry">

Chamberlain, Scott, Kara Woo, Andrew MacDonald, Naupaka Zimmerman, and
Gavin Simpson. 2025. *<span class="nocase">pangaear</span>: Client for
the ’Pangaea’ Database*. Manual.

</div>

<div id="ref-Diaconu.2020" class="csl-entry">

Diaconu, Andrei-Cosmin, Ioan Tanţău, Klaus-Holger Knorr, Werner Borken,
Angelica Feurdean, Andrei Panait, and Mariusz Gałka. 2020. “A
Multi-Proxy Analysis of Hydroclimate Trends in an Ombrotrophic Bog over
the Last Millennium in the Eastern Carpathians of Romania.”
*Palaeogeography, Palaeoclimatology, Palaeoecology* 538 (January):
109390. <https://doi.org/10.1016/j.palaeo.2019.109390>.

</div>

<div id="ref-Galka.2022a" class="csl-entry">

Gałka, Mariusz, Andrei-Cosmin Diaconu, Angelica Feurdean, Julie Loisel,
Henning Teickner, Tanja Broder, and Klaus-Holger Knorr. 2022. “Relations
of Fire, Palaeohydrology, Vegetation Succession, and Carbon
Accumulation, as Reconstructed from a Mountain Bog in the Harz Mountains
(Germany) During the Last 6200 Years.” *Geoderma* 424 (October): 115991.
<https://doi.org/10.1016/j.geoderma.2022.115991>.

</div>

<div id="ref-Galka.2022" class="csl-entry">

Gałka, Mariusz, Adam Hölzer, Angelica Feurdean, Julie Loisel, Henning
Teickner, Andrei-Cosmin Diaconu, Marta Szal, Tanja Broder, and
Klaus-Holger Knorr. 2022. “Insight into the Factors of Mountain Bog and
Forest Development in the Schwarzwald Mts.: Implications for Ecological
Restoration.” *Ecological Indicators* 140 (July): 109039.
<https://doi.org/10.1016/j.ecolind.2022.109039>.

</div>

<div id="ref-Reuter.2019b" class="csl-entry">

Reuter, Hendrik, Julia Gensel, Marcus Elvert, and Dominik Zak. 2019.
“FTIR, CuO Lignin, and Bulk Decomposition Data of a 75-Day Anoxic
Phragmites Australis Litter Decomposition Experiment in Soil Substrates
from Three Northeast German Wetlands.” PANGAEA.
<https://doi.org/10.1594/PANGAEA.902181>.

</div>

</div>
