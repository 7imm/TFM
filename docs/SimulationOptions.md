# Simulation Options

This document lists all options to choose models within the simulation framework
TFM. The models can be defined in a configuration file representing a Fortran
name list called `"tfm.conf"`. An example file can be found in
[`example/tfm.conf`](../example/tfm.conf).

---

## `solve_density`

#### `"false"`
If the argument `solve_density="false"` is passed, TFM does not solve for
densification.



#### `"herron1980"`
The classical densification model of Herron & Langway (1980).

**Herron, M. M. and Langway, C. C. (1980)**.
*Firn Densification: An Empirical Model.*
Journal of Glaciology 25.93, pp. 373-385.
DOI: https://doi.org/10.3189/S0022143000015239



#### `"li2003"`
The firn densification model of Li & Zwally (2003) is another modification of
the Herron & Langway (1980) model.

**Li, J., Zwally, H. J., Corneja, H., and Yi, D. (2003)**.
*Seasonal variation of snow surface elevation in North Greenland as modeled and
detected by satellite radar altimetry.*
Annals of Glaciology 37, pp. 223-238.
DOI: https://doi.org/10.3189/172756403781815889



#### `"helsen2008"`
The firn densification model of Helsen et al. (2008).

**Helsen, M. M., van den Broeke, M. R., ven den Wal, R. S. W.,
van de Berg, W. J., van Meijgaard, E., Davis, C. H., Li, Y.,
and Goodwin, I. (2008)**.
*Elevation Changes in Antarctica Mainly Determined by Accumulation Variability.*
Science 320.5883, pp. 1626-1629.
DOI: https://doi.org/10.1126/science.1153894



#### `"arthern2010"`
The firn densification model of Arthern et al. (2010) is especially tuned for
Antarctica. It is based on the approach of Herron & Langway (1980).

**Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R.
(2010)**.
*In situ measurements of Antarctic snow compaction compared with predictions of
models.*
Journal of Geophysical Research: Earth Surface 115.F3.
DOI: https://doi.org/10.1029/2009JF001306



#### `"ligtenberg2011"`
The firn densification model of Ligtenberg et al. (2011), following the approach
of Herron & Langway (1980).

**Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R. (2011)**.
*An improved semi-empirical model for the densification of Antarctic firn.*
The Cryosphere 5, pp. 809-819.
DOI: https://doi.org/10.5194/tc-5-809-2011



#### `"simsonsen2013`"
The firn densification model of Simonsen et al. (2013), following the approach
of Herron & Langway (1980).

**Simonsen, S. B. Stenseng, L., Adalgeisdottir, G., Fausto, R. S., Hvidberg, C.
S., and Lucas-Picher, P. (2013)**.
*Assessing a multilayered dynamic firn-compaction model for Greenland with
ASIRAS radar mesaurements.*
Journal of Glaciology 59.215, pp. 545-558.
DOI: https://doi.org/10.3189/2013JoG12J158



#### `"medley2020"`
The firn densification of Medley et al. (2020) is a tune of the model by
Herron & Langway (1980). Note that the implementation in TFM is based on the
2020 preprint of the paper. The final paper was published in the Cryosphere in
2022 and differs from the preprint in central aspects.

**Medley, B., Neumann, T. A., Zwally, H. J., and Smith, B. E. (2020)**.
*Forty-year Simulations of Firn Processes over the Greenland and Atnarctic Ice
Sheets.*
The Cryosphere Dicussions [preprint], in review.
DOI: https://doi.org/10.5194/tc-2020-266



#### `"arthern1998`"
The firn densification model of Arthern & Wingham (1998) establishes a process
based densification approach, based on sintering theory.

**Arthern, R. J. and Wingham, D. J. (1998)**.
*The Natural Fluctuations of Firn Densification and Their Effect on the Geodetic
Determination of Ice Sheet Mass Balance.*
Climate Change 40, pp. 605-624.
DOI: https://doi.org/10.1023/A:1005320713306



#### `"breant2017"`
Like the model of Arthern & Wingham (1998), the model of Breant et al. (2017)
incorporates aspects of sintering theory.

**Breant, C., Martinerie, P., Orsi, A., Arnaud, L., and Landais, A. (2017)**.
*Modelling firn thickness evolution during the last deglaciation: constraints on
sensitivity to temperature and impurities.*
Clim. Past 13, pp. 833-853.
DOI: https://doi.org/10.5194/cp-13-833-2017



#### `"sintering"`
The `"sintering"` model is an experimental process based model similar to that
of Arthern & Wingham (1998).



#### `"gagliardini1998"`
The model of Gagliardini & Meyssonnier (1998) is based on a so-called cell model
approach, which has its origin in theory of porous media and homogenization. It
is more complex than other models and requires iterative solving of a non-linear
constitutive equation.

**Gagliardini, O. and Meyssonnier, J. (1998)**.
*Flow simulation of a firn-covered cold glacier.*
Annals of Glaciology 24, pp. 242-248.
DOI: https://doi.org/10.3189/S0260305500012246



#### `"zwinger2007"`
The firn densificaiton model of Zwinger et al. (2007) is based on the model of
Gagliardini & Meyssonnier (1998). It differs from this model in certain aspects.

**Zwinger, T., Greve, R., Gagliardini, O., Shiraiwa, T.,
and Lyly, M. A. (2007)**.
*A full Stokes-flow thermo-mechanical model for firn and ice applied to the
Gorshkov crater glacier, Kamchatka.*
Annals of Glaciology 45, pp. 29-37.
DOI: https://doi.org/10.3189/172756407782282543



#### `"greve2009"`
In their textbook Greve & Blatter (2009) again describe the model of
Gagliardini & Meyssonnier (1998) with small modifications.

**Greve, R. and Blatter, H. (2009)**.
*Dynamics of Ice Sheets and Glaciers.*
Springer, Berlin.


#### `"timmsfit`"
Schultz (2024) extensively reviewed the model of Gagliardini & Meyssonnier
(1998) and modified it for a greater dataset of firn measurements. A description
of the model can be found in the thesis listed below.

**Schultz, T. (2024)**.
*"Physical Modeling of Firn: Densification, Temperature, Grain Growth, and Water
Retention"*.
Dissertation. TU Darmstadt. DOI:
https://doi.org/10.26083/tuprints-00027894.

---

## `solve_temperature`

#### `"false"`
If `"false"` is passed to `solve_temperature` temperature diffusion along the
firn profile is not computed.



#### `"true"`
If `solve_temperature="true"` temperature diffusion along the firn profile is
computed.

---

## `solve_heat_capacity`
TFM calculates the specific heat capacity of firn using a mixture theory
approach (see Schultz, 2024, "4.1.4 Specific Heat Capacity of Firn", p. 94 ff.).
The option `solve_heat_capacity` defines how the specific heat capacity of ice,
from which the specific heat capacity of firn is derived, is defined.



#### `"paterson1994"`
In his 1994 textbook, Paterson provides a constant value for the specific heat
capacity of ice.

**Paterson (1994)**.
*The Physics of Glaciers*.
Butterworth-Heinemann, Oxford, Amsterdam, Boston, London, New York, Paris,
San Diego, San Francisco, Singapore, Sydney, Tokyo. Third Edition.
ISBN 0 7506 4742 6



#### `"cuffey2010"`
In the 2010 edition of the book by Paterson, now Cuffey & Paterson (2010), a
temperature-dependent definition of the specific heat capacity of ice is
provided.

**Cuffey, K. M. and Paterson, W. S. B. (2010)**.
*The Physics of Glaciers*.
Butterworth-Heinemann, Amsterdam, Boston, Heidelberg, London, New York, Oxford,
Paris, San Diego, San Francisco, Singapore, Sydney, Tokyo. Fourth Edition.
ISBN 978-0-12-369461-4

---

## `solve_thermal_conductivity`
The option `solve_thermal_conductivity` defines how the thermal conductivity of
dry firn is computed. For the computation of the thermal conductivity of wet
firn, see `solve_liquid_thermal_conductivity`.



#### `"sturm1997"`
The parameterization by Sturm et al. (1997) calculates the thermal conductivity
of snow based on the density.

**Sturm, M., Holmgren, J., König, M., and Morris, K. (1997)**.
**The thermal conductivity of seasonal snow.**
Journal of Glaciology 43.143, pp. 26-41.
DOI: https://doi.org/10.3189/S0022143000002781



#### `"marchenko2019"`
The parameterization of the thermal conductivity by Marchenko et al. (2019)
relies solely on the density of the material. 

**Marchenko, S., Cheng, G., Löstedt, P., Pohjola, V., Pettersson, R.,
van Pelt, W., and Reijmer, C. (2019)**.
*Thermal conductivity of firn at Lomonosovfonna, Svalbard, derived from
subsurface temperature measurements.*
The Cryosphere 13.7, pp. 1843-1859.
DOI: https://doi.org/10.5194/tc-13-1843-2019



#### `"calonne2019"`
The description of the thermal conductivity by Calonne et al. (2019).

**Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F.,
and Geindreau, C. (2019)**.
*Thermal Conductivity of Snow, Firn, and Porous Ice From 3-D Image-Based
Computations.*
Geophysical Research Letters 46.22, pp. 13079-13089.
DOI: https://doi.org/10.1029/2019GL085228



#### `"miller1969upperbound"`
The model of Miller (1969) describes the thermal conductivity of a porous
material using a mixture theory approach, based on the material's porosity and
the thermal conductivity of its phases. The option `"miller1969upperbound"`
computes the upper bound provided by Miller (1969).

**Miller, M. N. (1969)**.
*Bounds for Effective Electrical, Thermal, and Magnetic Properties of
Heterogeneous Materials.*
Journal of Mathematical Physics 10.11, pp. 1988-2004.
DOI: https://doi.org/10.1063/1.1664794



#### `"miller1969lowerbound"`
The model of Miller (1969) describes the thermal conductivity of a porous
material using a mixture theory approach, based on the material's porosity and
the thermal conductivity of its phases. The option `"miller1969lowerbound"`
computes the lower bound provided by Miller (1969).

**Miller, M. N. (1969)**.
*Bounds for Effective Electrical, Thermal, and Magnetic Properties of
Heterogeneous Materials.*
Journal of Mathematical Physics 10.11, pp. 1988-2004.
DOI: https://doi.org/10.1063/1.1664794



#### `"geometricmean"`
The option `"geometricmean"` applies volume based geometric mean weighting to
the different phases of the porous material firn to compute the effective
thermal conductivity.

---

## `solve_liquid_thermal_conductivity`

The computation of the thermal conductivity of wet firn is based on the
suggestions made by Kaviany (1991, p. 126) in his textbook. He suggests to use
either geometric mean weighting or a Voigt weighting model to compute the
effective thermal conductivity of a three phase medium like wet firn, consisting
of ice, air, and water.
To do this, the thermal conductivity at water saturation has to be known, see
`solve_saturation_thermal_conductivity`.

**Kaviany, M. (1991)**.
*Principles of Heat Transfer in Porous Media.*
in: Ling, F. F. (edt.) Mechanical Engineering Series, Springer-Verlag, New York,
Berlin, Heidelberg.
DOI: https://doi.org/10.1063/1.1664794



#### `"false"`
If the option `"false"` is passed to `solve_liquid_thermal_conductivity` the
effective thermal conductivity of wet firn is not calculated. The thermal
conductivity of dry firn is used instead.


#### `"geometricmean"`
Volumetric geometric mean weighting to compute the effective thermal
conductivity of wet firn, based on the thermal conductivity of dry firn and the
thermal conductivity of firn at water saturation.


#### `"voigt"`
Voigt model weighting to compute the effective thermal conductivity of wet firn,
based on the thermal conductivity of dry firn and the thermal conductivity of
firn at water saturation.

---

## `solve_saturation_thermal_conductivity`

The thermal conductivity of water saturated firn is needed to compute the
effective thermal conductivity of firn at lower water saturation levels using
mixture theory. Therefore, if `solve_liquid_thermal_conductivity` is given,
`solve_saturation_thermal_conductivity` must not be `"false"`. There is no
certain method for computing the effective thermal conductivity of the two
phase material consisting of ice and water. However, various mixture theory
weighting approach are provided by TFM.



#### `"false"`
The option `"false"` can only be passed to
`solve_saturation_thermal_conductivity` if `solve_liquid_thermal_conductivity`
is also `"false"`.



#### `"geometricmean"`
Geometric mean weighting to compute the effective thermal conductivity of water
saturated firn, based on the thermal conductivities of ice and water at the
pressure melting point.


#### `"voigt"`
Voigt model weighting to compute the effective thermal conductivity of water
saturated firn, based on the thermal conductivities of ice and water at the
pressure melting point.

**Voigt, W. (1889)**.
*Ueber die Beziehung zwsichen den beiden Elasticitätsconstanten Isotroper
Körper.*
Annalen der Physik 274.12, pp. 573-587.
DOI: https://doi.org/10.1002/andp.18892741206



#### `"reuss"`
Reuss model weighting to compute the effective thermal conductivity of water
saturated firn, based on the thermal conductivities of ice and water at the
pressure melting point.

**Reuss, A. (1929)**.
*Berechnung der Fließgrenze von Mischkristallen auf Grund der
Plastizitätbedingung für Einkristalle.*
ZAMM - Journal of Applied Mathematics and Mechanics 9.1, pp. 49-58.
DOI: https://doi.org/10.1002/zamm.19290090104



#### `"miller1969upperbound"`
Upper bound method provided by Miller (1969) to compute the effective thermal
conductivity of water saturated firn, based on the thermal conductivities of ice
and water at the pressure melting point.

**Miller, M. N. (1969)**.
*Bounds for Effective Electrical, Thermal, and Magnetic Properties of
Heterogeneous Materials.*
Journal of Mathematical Physics 10.11, pp. 1988-2004.
DOI: https://doi.org/10.1063/1.1664794



#### `"miller1969lowerbound"`
Lower bound method provided by Miller (1969) to compute the effective thermal
conductivity of water saturated firn, based on the thermal conductivities of ice
and water at the pressure melting point.

**Miller, M. N. (1969)**.
*Bounds for Effective Electrical, Thermal, and Magnetic Properties of
Heterogeneous Materials.*
Journal of Mathematical Physics 10.11, pp. 1988-2004.
DOI: https://doi.org/10.1063/1.1664794

---

## `solve_liquid`



#### `"false"`
If `"false"` is passed to `solve_liquid`, water flow through the firn column is
not solved. Water infiltration is ignored.



#### `"richardsequation"`
Option to solve Richards Equation for unsaturated water flow through the firn
column. If the option `"richardsequation"` is passed to `solve_liuqid`, the
option `solve_van_genuchten` cannot be `"false"`.



#### `"bucket"`
With the option `"bucket"` passed to `solve_liquid` a simple tipping bucket
scheme is used to simulate water flow through the firn column.

NOTE: This method has not been tested in detail. Use it only with great caution
and not in production.

---

## `solve_van_genuchten`



#### `"daanen2009"`
The parameterization of the van Genuchten parameters for snowpack by
Daanen & Nieber (2009) relies solely on the grain radius.

**Daanen, R. P. and Nieber, J. L. (2009)**.
*Model for Coupled Liquid Water Flow and Heat Transport with Phase Change in a
Snowpack.*
Journal of Cold Regions Engineering 23.2, pp. 43-68.
DOI: https://doi.org/10.1061/(ASCE)0887-381X(2009)23:2(43)



#### `"yamaguchi2010"`
The Yamaguchi et al. (2010) paraterization of the water retention curve for snow
is based, like that of Daanen & Nieber (2009), on the material's grain radius.

**Yamaguchi, S., Katsushima, T., Sato, A., and Kumakura, T. (2010)**.
*Water retention curve of snow with different grain sizes.*
Cold Regions Science and Technology 64.2, pp. 87-93.
DOI: ttps://doi.org/10.1016/j.coldregions.2010.05.008



#### `"yamaguchi2012"`
In contrast to the parameterization of Daanen & Nieber (2009) and Yamaguchi
et al. (2010) the paramterization of the van Genuchten parameters by Yamaguchi
et al. (2012) considers the density of snow as well as its grain size.

**Yamaguchi, S., Watanae, K., Katsushima, T., Sato, A.,
and Kumakura, T. (2012)**.
*Dependence of the water retention curve of snow on snow characteristics.*
Annals of Glaciolog 51.61, pp. 6-12.
DOI: https://doi.org/10.3189/2012AoG61A001

---

## `solve_grain_growth`



#### `"false"`
If `"false"` is passed to `solve_grain_growth`, grain growth is not solved
within the simulation. The grain radius will be defined by the initial and
boundary conditions.



#### `"arthern2010"`
The grain growth model presented by Arthern et al. (2010) is intended for use in
Antarctic dry conditions.

**Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R.,
and Thomas, E. R. (2010)**.
*In situ measurements of Antarctic snow compaction compared with predictions of
models.*
Journal of Geophysical Research: Earth Surface 115.F3.
DOI: https://doi.org/10.1029/2009JF001306



#### `"zwally2002"`
Grain growth follwing Zwally & Li (2002) as described in Li & Zwally (2011).

**Zwally, H. J. and Li, J. (2002)**.
*Seasonal and interannual variations of firn densification and ice-sheet
elevation at the Greenland summit.*
Journal of Glaciology 48.171.
DOI: https://doi.org/10.3189/172756502781831403

**Li, J. and Zwally, H. J. (2011)**.
*Modeling of firn compaction for estimating ice-sheet mass change from observed
ice-sheet elevation change.*
Annals of Glaciology 52.59.
DOI: https://doi.org/10.3189/172756411799096321



#### `"brun1989"`
The grain growth model by Brun (1989) describes grain growth with respect to the
water content within snow.

**Brun, E. (1989)**.
*Investigation On Wet-Snow Metamorphism in Respect of Liquid-Water Content.*
Annals of Glaciology 13, pp. 22-26.
DOI: https://doi.org/10.3189/S0260305500007576



#### `"tusima1978"`
The grain growth model by Tusima (1978), like the one by Brun (1989), considers
enhanced grain growth in aquatic conditions. The original paper is wirtten in
Japanese. However, the paper by Katsushima et al. (2009) describes the model.

Grin growth following Tusima (1978) as described by Katsushima et al. (2009)
(the original paper is in Japanese).

**Tusima, K. (1978)**.
*Grain coarsening of ice particles immersed in pure water.*
Sppyo 40.4, pp. 155-165.

**Katsushima, T., Kamakura, T., and Takeuchi, Y. (2009)**.
*A multiple layer model including a parameterization of vertical water channel
process in snowpack.*
Cold Regions Science and Technology 59, pp. 143-151.
DOI: https://doi.org/10.1016/j.coldregions.2009.09.002



#### `"katsushima2009"`
Grain growth following Katsushima et al. (2009) combining the parameterization
of Brun (1989) and Tusima (1978) depending on the liquid water content of firn.

**Katsushima, T., Kamakura, T., and Takeuchi, Y. (2009)**.
*A multiple layer model including a parameterization of vertical water channel
process in snowpack.*
Cold Regions Science and Technology 59, pp. 143-151.
DOI: https://doi.org/10.1016/j.coldregions.2009.09.002