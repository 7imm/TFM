```
   _\/  \/_      ______________________  ___
    _\/\/_       ___  __/__  ____/__   |/  /
_\_\_\/\/_/_/_   __  /  __  /_   __  /|_/ /
 / /_/\/\_\ \    _  /   _  __/   _  /  / /
    _/\/\_       /_/    /_/      /_/  /_/
    /\  /\
```

TFM is the firn simulation framework that I developed as part of my PhD thesis.
It is a classical 1D firn model that can simulate:

- firn densification
- temperature evolution (including specific heat capacity and thermal
conductivity in the presence of liquid water)
- unsaturated water flow through the firn body by solving Richards Equation
- water retention due to refreezing
- grain growth

The key concepts of the model are described in the thesis:

**Schultz, T (2024)**.
*"Physical Modeling of Firn: Densification, Temperature, Grain Growth, and Water
Retention"*.
Dissertation. TU Darmstadt. DOI:
https://doi.org/10.26083/tuprints-00027894.

If you want to use TFM and have a question, simply write me a
[mail](mailto:tfm@firn.cool).



# Quickstart
TFM is written in modern Fortran. For convenience, I decided to release it as a
[fpm](https://fpm.fortran-lang.org/index.html) package, making it easy to use
even for non-Fortran programmers. This also allows including TFM easily in own
projects. To set up the Fortran package manager, see the
[official documentation](https://fpm.fortran-lang.org/install/index.html).



## Installation
Clone [this GitHub repository](https://github.com/7imm/TFM) and change into the
cloned directory.

```
git clone https://github.com/7imm/TFM
cd TFM
```

To install TFM run

```
fpm install --flag -fbackslash
```

The option `--flag -fbackslash` is not necessary but ensures that the progress
indicator is working correctly.

To run the generic example case coming with the repository change into the
`example` directory, create an output directory called `tfm_output`, and run the
model calling `TFM`.

```
cd example
mkdir tfm_output
TFM
```

Congratulations! You have done your first simulation using TFM.



# How to use TFM
In general, TFM is a simulation framework that can be used in many ways. It is
meant to be used within your own project. However, the repository comes with an
example you can find in `app/main.f90`. This program can serve as a starting
point for your own project. All routines used by `app/main.f90` can be found in
`src/`.

Initially, all things needed to run TFM are a project file like `app/main.f90`,
a configuration file defining what the model shall do, a file defining the
initial firn profile, and an input file prescribing the upper boundary
conditions.

The example files coming with this repository are meant to be as simple as
possible. No fiddling around with `netcdf`, `json` or `xml` files, but simple
`txt` files. However, you can, of course, write your own scripts for reading
such files and passing their data to TFM. The next sections will describe the
input files for the example case.

I want to illustrate how to use TFM using the example you can find in the
directory `example`. There you will find three files:

- [`example/tfm.conf`](example/tfm.conf)
- [`example/example_init.dat`](example/example_init.dat)
- [`example/example_forcing.dat`](example/example_forcing.dat)

In the following, I will explain these files. You can find a more detailed
description of the example case in [`docs/Example.md`](docs/Example.md).



## The configuration file
The configuration file for TFM is  hardcoded and called `tfm.conf`. If you want
to change this name, you can find it in [`app/main.f90`](app/main.f90) or you
can create your own project with another name. For now, we stick with the
example file.
```
&CONFIG
  solve_density="arthern2010",
  solve_temperature="true",
  solve_heat_capacity="cuffey2010",
  solve_thermal_conductivity="calonne2019",
  solve_liquid_thermal_conductivity="geometricmean",
  solve_saturation_thermal_conductivity="geometricmean",
  solve_liquid="richards_equation",
  solve_van_genuchten="yamaguchi2012",
  solve_grain_growth="katsushima2009",
  forcing_input_file="./example_forcing.dat",
  initial_input_file="./example_init.dat",
  time_step=86400.0,
  spinup=-1,
  max_profile_length=-1.0,
  max_profile_age=-1.0,
 /
```

The configuration file contains a Fortran name list. Such a name list can be
parsed in Fortran easily. However, the order of the options cannot be changed.
It defines the different models that shall be used in the simulation, as well as
some other things needed for the simulation, like the time step.

For a complete list of all options, see
[`docs/SimulationOptions.md`](docs/SimulationOptions.md).

---

#### `solve_density`
This option defines which densification model shall be used during the
simulation. In the example case, the model of Arthern et al. (2010) is used. The
option has to be a string. The string `"false"` is valid. In this case,
densification will not be solved. If a non-valid string is given, the model will
throw an error and will stop working. For a full list of possible options and
the corresponding references, see
[`docs/SimulationOptions.md`](docs/SimulationOptions.md).

**Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R.
(2010)**.
*In situ measurements of Antarctic snow compaction compared with predictions of
models.*
Journal of Geophysical Research: Earth Surface 115.F3.
DOI: https://doi.org/10.1029/2009JF001306

---

#### `solve_temperature`
The options for `solve_temperature` are either `"true"` or `"false"`. Either
temperature diffusion along the firn profile will be solved or not.

---

#### `solve_heat_capacity`
In TFM, the specific heat capacity of firn is calculated using a mixture theory
approach. This is especially relevant if liquid water in the firn profile
refreezes (see Schultz, 2024, pp. 94 ff.). The option `solve_heat_capacity`
determines how the specific heat capacity of ice is determined. There are only
two options, which are `"paterson1994"` and "`cuffey2010`". While Paterson
(1994) gave a constant value in his 1994 textbook, the later edition by
Cuffey & Paterson (2010) provides a temperature-dependent relation.

**Paterson (1994)**.
*The Physics of Glaciers*.
Butterworth-Heinemann, Oxford, Amsterdam, Boston, London, New York, Paris,
San Diego, San Francisco, Singapore, Sydney, Tokyo. Third Edition.
ISBN 0 7506 4742 6

**Cuffey, K. M. and Paterson, W. S. B. (2010)**.
*The Physics of Glaciers*.
Butterworth-Heinemann, Amsterdam, Boston, Heidelberg, London, New York, Oxford,
Paris, San Diego, San Francisco, Singapore, Sydney, Tokyo. Fourth Edition.
ISBN 978-0-12-369461-4

---

#### `solve_thermal_conductivity`
The option `solve_thermal_conductivity` defines how the thermal conductivity of
dry firn is computed during the simulation. In the present case, the
parameterization of Calonne et al. (2019) is used. There are many such
parameterizations for snow and firn. Only a few of them are implemented in TFM.
Another option to compute the thermal conductivity of firn is again to use a
mixture approach, based on the thermal conductivities of ice and air. For a
complete list of options, see
[`docs/SimulationOptions.md`](docs/SimulationOptions.md).

**Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F.,
and Geindreau, C. (2019)**.
*Thermal Conductivity of Snow, Firn, and Porous Ice From 3-D Image-Based
Computations.*
Geophysical Research Letters 46.22, pp. 13079-13089.
DOI: https://doi.org/10.1029/2019GL085228

---

#### `solve_liquid_thermal_conductivity`
While there are many parametrizations for the thermal conductivity of dry firn,
I'm not aware of any models for wet firn describing the thermal conductivity.
Kaviany (1991, p. 126) suggests in his textbook to use either geometric mean
weighting or a Voigt model. Hence, `"geometricmean"` and `"voigt"` are the two
options that can be used to calculate the thermal conductivity of wet firn. A
third option is to pass `"false"`. In this case, the mixture theory approach is
ignored, and the thermal conductivity of dry firn is used for the simulation.

**Kaviany, M. (1991)**.
*Principles of Heat Transfer in Porous Media.*
in: Ling, F. F. (edt.) Mechanical Engineering Series, Springer-Verlag, New York,
Berlin, Heidelberg.
DOI: https://doi.org/10.1063/1.1664794

---

#### `solve_saturation_thermal_conductivity`
To solve for the thermal conductivity using a mixture theory approach, the
thermal conductivity of firn at water saturation is needed. If 
`solve_liquid_thermal_conductivity` is defined and
`solve_saturation_thermal_conductivity="false"` TFM throws an error. The option
allows defining a method to compute this property. This is also done using
mixture theory. In the given case, geometric mean weighting is also applied. For
a detailed description of how the mixture theory approach in TFM works, it is
referred to the corresponding section in Schultz (2024, see Sec. 4.1.2
"Thermal Conductivity of Firn in the Presence of Liquid Water", p. 90 ff.). For
a full list of options for `solve_saturation_thermal_conductivity` see
[`docs/SimulationOptions.md`](docs/SimulationOptions.md).

---

#### `solve_liquid`
The option `solve_liquid` defines if and how water flow through the firn column
is solved. If `solve_liquid="false"` water transport is not simulated. The other
options are `solve_liquid="richardsequation"`, like in the example case, and
`solve_liquid="bucket"`. Depending on the given option, either Richards Equation
is solved to simulate unsaturated flow or a bucket scheme is employed. Note, the
bucket scheme was never really tested in detail. It is not recommended to use
this option in production! If the Richards Equation shall be solved, the option
`solve_van_genuchten` has to be given.

---

#### `solve_van_genuchten`
The option `solved_van_genuchten` determines which parameterization of the van
Genuchten parameters for snow/firn shall be used to model the water retention
curve of the porous material. In the given case, the parameterization of
Yamaguchi et al. (2012) is used. Other options are the models of Yamaguchi
(2010) and Daanen & Nieber (2009). For full references, see
[`docs/SimulationOptions.md`](docs/SimulationOptions.md).

**Yamaguchi, S., Watanae, K., Katsushima, T., Sato, A.,
and Kumakura, T. (2012)**.
*Dependence of the water retention curve of snow on snow characteristics.*
Annals of Glaciology 51.61, pp. 6-12.
DOI: https://doi.org/10.3189/2012AoG61A001

---

#### `solve_grain_growth`
To define the grain growth model that shall be solved during the simulation, the
option `solve_grain_growth` has to be defined. If `solve_grain_growth="false"`
grain growth is not solved. Different models either incorporating abnormal grain
growth due to liquid water or not are available. For a complete list including
corresponding references, see
[`docs/SimulationOptions.md`](docs/SimulationOptions.md). In the example case,
the model defined by Katsushima et al. (2009) is employed.

**Katsushima, T., Kamakura, T., and Takeuchi, Y. (2009)**.
*A multiple layer model including a parameterization of vertical water channel
process in snowpack.*
Cold Regions Science and Technology 59, pp. 143-151.
DOI: https://doi.org/10.1016/j.coldregions.2009.09.002

---

#### `forcing_input_file`
The option `forcing_input_file` is a string defining where to find the forcing
file for the simulation. In the example case, the forcing file is called
`"example_forcing.dat"`. It can be found in
[`example/example_forcing.dat`](example/example_forcing.dat). For information
about how the file has to be formatted, see
["The forcing file"](## The forcing file).

---

#### `initial_input_file`
Like `forcing_input_file`, `inital_input_file` defines where to find the file
that defines the initial conditions of the simulated firn profile. In the
example case, the initial input file is found at
[`example/example_init.dat`](example/example_init.dat). For information about
formatting the init file, see ["The init file"](## The init file)

---

#### `time_step`
The option `time_step` defines the time step of the simulation. It is given as a
floating-point number, or `REAL` in Fortran. It defines the time step in seconds.
In the example case, the time step is `time_step=86400.0` which is the number of
seconds in a day.

---

#### `spinup`
The option `spinup` allows doing a spin-up simulation before the actual
simulation. The option is given as a signed integer, defining how many years the
spin-up should last. If `spinup=0` or the given integer is negative, like in the
example case, no spin-up is computed.

Let's assume `spinup=100`. Then a hundred-year spin-up is computed using the
given time step. The initial conditions of the firn profile are defined by
`initial_input_file`. The forcing for the spin-up is computed as a means of the
forcing provided via `forcing_input_file`. After `365 * 100 = 36500` times
steps, the actual simulation would start using the forcing provided via
`forcing_input_file`.

Using a spin-up is of particular interest if no well-defined initial conditions
are available. For example, one could define a ten-meter-long initial profile
consisting entirely of ice by defining an initial density of
`density = ICE_DENSITY = 917.0`. A spin-up would then result in a firn profile
representing mean climatic conditions from the forcing file.

---

#### `max_profile_length`
TFM follows a Lagrangian simulation approach, where layers are added at the top
of the firn profile if positive solid accumulation is present within a time
step. This means that the number of layers within the firn profile is usually
growing. This in turn leads to a slowdown of the model during the simulation.

As TFM is fast, I decided not to merge layers, like it is done in other firn
simulation frameworks. TFM always works at full resolution. However, often deep
layers of the firn column, near ice density, are not of so much interest.
Therefore, the option `max_profile_length` allows cutting the firn profile at
the bottom. `max_profile_length` accepts a signed floating-point number. If this
number is negative, like in the example case, the profile is not cut at the
bottom. If we assume `max_profile_length=100.0` and the firn profile exceeds a
total length of 100 meters, the corresponding layers are deleted from the bottom
of the profile.

---

#### `max_profile_age`
Like `max_profile_length`, `max_profile_age` allows limiting the length of a
simulated firn profile by its maximum age. Again, the option accepts a signed
floating-point number representing the maximum age in years. If this number is
negative, the option is ignored.



## The init file
Next, we take a look at the init file. The first five lines and the last five
lines of [`example/example_init.dat`](example/example_init.dat) are

```
-50.00 914.89 261.15 0.00463 0.00 30214099834
-49.90 914.86 261.15 0.00462 0.00 30141826367
-49.80 914.83 261.15 0.00461 0.00 30069555232
-49.70 914.80 261.15 0.00460 0.00 29997286463
-49.60 914.77 261.15 0.00459 0.00 29925020091
[...]
-0.40 368.24 261.15 0.00052 0.00 115143730
-0.30 366.17 261.15 0.00051 0.00 86135133
-0.20 364.11 261.15 0.00051 0.00 57289587
-0.10 362.05 261.15 0.00050 0.00 28606782
-0.00 360.00 261.15 0.00050 0.00 86400
```

It is a simple `txt` file representing a table with the following columns.

|      | depth $z$    | density $\rho$          | temperature $T$ | grain radius $r$ | liquid water content $\theta$ | age $\chi$   |
|------|--------------|-------------------------|-----------------|------------------|-------------------------------|--------------|
| unit | $\mathrm{m}$ | $\mathrm{kg \, m^{-3}}$ | $\mathrm{K}$    | $\mathrm{m}$     | $\mathrm{1}$                  | $\mathrm{s}$ |

The properties are self-explanatory. TFM usually uses SI units (with some
exceptions like `max_profile_age`). This sometimes leads to great numbers, like
in the case of age, or small numbers, like in the case of accumulation rate.
However, I decided to stick to SI units because I'm not a fan of derived units.

The depth column accepts negative and position floating-point numbers. For
example, a profile can start at ${z_\mathrm{bottom} = -25 \, \mathrm{m}}$ and
end at ${z_\mathrm{top} = +25 \, \mathrm{m}}$. TFM, unlike other models, does
not shift the top or bottom of the profile to a certain reference point.
However, the reference point of the profile is always the bottom layer.

Note, that the order of the file ranges from bottom to top. This is unusual when
it comes to data of firn profiles, but necessary because of how TFM works.



## The forcing file
Like the init file, the forcing file is a simple `txt` representing  a table. In
the following, you can see the first five lines and the last five lines of the
example file [`example/example_forcing.dat`](example/example_forcing.dat).

```
0 0 261.15 360.00 2.853881e-09 0.000000e+00 0.00050
86400 86400 261.41 360.00 2.853646e-09 0.000000e+00 0.00050
172800 172800 261.67 360.00 2.852942e-09 0.000000e+00 0.00050
259200 259200 261.92 360.00 2.851768e-09 0.000000e+00 0.00050
345600 345600 262.18 360.00 2.850124e-09 0.000000e+00 0.00050
[...]
314928000 314928000 259.86 360.00 2.848012e-09 0.000000e+00 0.00050
315014400 315014400 260.12 360.00 2.850124e-09 0.000000e+00 0.00050
315100800 315100800 260.38 360.00 2.851768e-09 0.000000e+00 0.00050
315187200 315187200 260.63 360.00 2.852942e-09 0.000000e+00 0.00050
315273600 315273600 260.89 360.00 2.853646e-09 0.000000e+00 0.00050
```

The columns store the following information:

|      | world time $t_w$ | simulation time $t$ | surf. temperature $\mathrm{T_0}$ | surf. density $\rho_0$  | solid accumulation $\dot{b}_s$ | liquid accumulation $\dot{b}_l$ | surf. grain radius $r_0$ |
|------|------------------|---------------------|----------------------------------|-------------------------|--------------------------------|---------------------------------|--------------------------|
| unit | $\mathrm{s}$     | $\mathrm{s}$        | $\mathrm{K}$                     | $\mathrm{kg \, m^{-3}}$ | $\mathrm{m \, weq. \, s^{-1}}$ | $\mathrm{m \, weq. \, s^{-1}}$  | $\mathrm{m}$             |

Each line of the table holds the information on one time step. TFM reads the
forcing file and interprets the number of lines in the file as the number of
time steps. Together with the time step, defined in the
[configuration file](## The configuration file) the total simulation time
follows.

The first two columns, "world time" and "simulation time", are not used
effectively by TFM. I implemented them because I thought it would be useful for
certain simulations to match the real-world time to the simulation time.
However, with the current [output routine](## Output) world time and simulation
time are not written to output.

Surf. temperature, surf. density, and surf. grain radius are the prescribed
temperature, density, and grain radius at the firn profile surface for each time
step, respectively.

I decided to implement the rather abstract concepts of solid and liquid
accumulation because determining how much solid and liquid accumulation is
present over the simulation time is usually a preprocessing step. While other
models implement concepts like a degree-day model and solve these concepts at
runtime, I believe that this introduces unnecessary complexity to a firn model.
The same is true for other properties, like the surface grain radius, that may
depend on the climatic conditions at the location of the firn profile.

The solid and liquid accumulation rates are given in units of meter water
equivalent per second. This unit is quite unusual, but sticks to the paradigm to
use SI-units. However, it is of course a compromise. What do solid accumulation
and liquid accumulation mean? Solid accumulation can be positive and negative.
Positive solid accumulation means that snow is added at the top of the firn
profile. It represents common accumulation due to snowfall. If the solid
accumulation becomes negative, mass is removed from the top of the firn profile.
This can be due to several reasons, of which the most common one is probably
surface melting. Surface melting, of course, leads to the production of melt
water, which brings us to liquid accumulation. Liquid accumulation is always
positive. It describes how much liquid water is infiltrated into the firn
column. This can be due to meltwater production at the surface, but, for
example, also due to rainfall.



## Output
For example, case TFM writes an output file for every time step in a directory
called `tfm_output`. If this directory does not exist, TFM will crash. The name
of the directory is hard-coded in [`app/main.f90`](app/main.f90). The naming of
the output files follows the scheme `tfm%06i.out`. Every output file contains a
table with the following firn profile properties.

|      | depth $z$    | density $\rho$          | temperature $T$ | grain radius $r$ | age $\chi$   | liquid water content $\theta$ |
|------|--------------|-------------------------|-----------------|------------------|--------------|-------------------------------|
| unit | $\mathrm{m}$ | $\mathrm{kg \, m^{-3}}$ | $\mathrm{K}$    | $\mathrm{m}$     | $\mathrm{s}$ | $1$                           |

Every line in the file represents a layer of the firn profile. Note that the
number of layers varies from time step to time step. The length of the output
file also varies.



# Overview: source files
You can find the source files of TFM in the `src/` directory. The following
provides a short overview of these files and what to find in them.

#### `src/tfm_constants.f90`
This module holds various common physical constants needed for the simulations.
However, many special constants are only needed once and are hard-coded in the
routines where they are needed.

#### `src/tfm_io.f90`
`tfm_io.f90` contains the module `tfm_io`, obviously including some routines for
IO. However, IO in TFM is not very elaborated and kept very simple.

#### `src/tfm_llstructure.f90`
Fortan has static and strong typing, which is one reason it's so fast. However,
TFM uses a Lagrangian model approach where layers are added to or deleted from
the simulated firn profile at every time step. Therefore, the amount of needed
memory is not always known at compile time. To be a little more flexible when it
comes to memory allocation, TFM uses a linked list structure to save the
profile's properties. This structure and related routines can be found in the
module `tfm_llstructure`.

#### `src/tfm_temperature.f90`
`tfm_temperature.f90` contains everything related to the simulation of
temperature in firn (in fact, almost everything, temperature change due to
refreezing is handled in `tfm_liquid`). This includes temperature diffusion but
also the computation of heat capacity and thermal conductivity for dry firn and
for firn in the presence of liquid water.

#### `src/tfm_density.f90`
The various densification models for firn can be found in `tfm_density.f90`. The
file holds more than one module, but all are related to densification.

#### `src/tfm_grain.f90`
Grain growth models are implemented in `tfm_grain.f90`. Besides the classical
grain growth models, the module implements routines for explicit and implicit
calculation of grain radius change from change of cross-sectional area and
volume.

#### `src/tfm_liquid.f90`
`tfm_liquid.f90` implements the scheme for solving Richards Equation and the
related van Genuchten model. This includes different parameterizations of van
Genuchten parameters for snow/firn. Numerical aspects like adaptive time
stepping for solving Richards Equation can also be found like a bucket scheme.
However, this bucket scheme is not tested very well and should be used with
great caution.

#### `src/tfm_num.f90`
`tfm_num.f90` is the heart of TFM. It's where all models are linked together in
a staggered scheme and boundary conditions are applied.

#### `src/tfm_essentials.f90`
This module holds some routines that fit nowhere else.



# How to Cite
If you use TFM, please cite my dissertation.

**Schultz, T (2024)**.
*"Physical Modeling of Firn: Densification, Temperature, Grain Growth, and Water
Retention*".
Dissertation. TU Darmstadt.
DOI: https://doi.org/10.26083/tuprints-00027894.


# License
Copyright 2024 Timm Schultz

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
