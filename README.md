# Oil palm growth and yield model in Julia

## Overview

Sawit.jl is a semi-mechanistic oil palm growth and yield model. Sawit.jl is a significant update to its predecessor [PySawit](https://github.com/cbsteh/PySawit), our early Python-based model. Because Sawit.jl is written in Julia and its code heavily optimized, this model runs very much faster by as much as 30 times than PySawit.

Sawit.jl was developed for Malaysian growing conditions. Moreover, its simulations were validated against historical data of field measurements from seven oil palm plantations in the country. These seven sites, across Peninsular Malaysia, differed in their soil type (Inceptisols and Ultisols), planting density (82 to 299 palms/ha), soil texture (27 to 59% clay and 7 to 67% sand), and annual rainfall amount (1800 to 2800 mm/yr).

The model showed overall good accuracy in simulating oil palm parameters (except for trunk weight) across diverse conditions, with model agreement metrics ranging from 6 to 27% for model absolute errors, -22 to +17% for model bias, and 0.38 to 0.98 for the Kling-Gupta Efficiency (KGE) index.

The model further predicted oil palm yield response to sudden rainfall changes, like during El Niño and La Niña events, and captured the effects of soil texture, rainfall, and meteorological factors on water deficits and crop photosynthesis.

## Usage
Clone this repo by

```
https://github.com/cbsteh/Sawit.jl.git
```

The source code is in the `src` folder, where the `core` subfolder comprises the model's six core components: (1) main loop, (2) meteorology, (3) photosynthesis, (4) energy balance, (5) soil water, and (6) crop growth.

The plotting functions are in the `plots` subfolder, and the utility functions in `util` subfolder.

Refer to `main.jl` for a quick way to run the model for a given site.

## Run the model
You can run the model in several ways, one of which is as follows:

```
using Sawit

jsonfname = "data/merlimau/input.json"
res = start(jsonfname)
plot_annual(res.annl, res.mi; saveplot=false)
plot_daily(res.daly, res.mi; saveplot=false)
```

where `jsonfname` holds the file name of the model input for the Merlimau site. This file is located in `data/merlimau` folder, and the name of the model input file is `input.json`. This file, in JSON plain text format, lists all the model parameters the model requires.

`res = start(jsonfname)` starts the simulation and upon completion, `res` contains the model daily and annual results, which can be retrieved by `res.daly` and `res.annl`, respectively. `res.mi` holds the model inputs, as read in earlier from the model input file.

`plot_annual` and `plot_daily` are functions to plot the model results on an annual and daily basis, respectively. Argument `saveplot` decides whether to save the plotted charts as a `PNG` picture format.

You can save the model output results to a `CSV` file as follows:

```
using CSV

CSV.write(“annual-results.csv”, res.annl)       # annual results to file
CSV.write(“daily-results.csv”, res.daly)           # daily results to file
```

which is possible because `res.annl` and `res.daly` return a `DataFrame` type.

## Input file (JSON format)
Model inputs are stored in a plain text file (JSON format) that can be easily edited in any text editor software. The model input file has the following parameters:

- `random generator seed`: seed for the random number generator (enter -1 for stochastic runs, else any integer value >0 for deterministic runs)
- `no. of simulation days`: number of simulation days (days), e.g., for 10 years, enter "3650"
- `data folder`: model input JSON file path. If code "@@" is given, the default working folder path will be used. This default folder path is the location of the model input `JSON` file. It is assumed the weather file and boxcar file (if any) will be located in this file path. The model output results will also be saved here.
- `boxcar filename`: name of Boxcar file, typically left blank
- `weather filename`: name of daily weather file
- `output filename`: name of file where model output will be saved. "-daily" and "-annual" will be appended to the given file name for daily and annual output, respectively
- `co2_path`: see comments for the `CCfn` function in the `mi.cc.jl` file, usually set as "past" for historical ambient CO2 levels
- `ta_path`: see comments for the `CCfn` function in the `mi.cc.jl` file. If set to "0", air temperature in the weather file will be read without change.
- `wind_path`: see comments for the `CCfn` function in the `mi.cc.jl` file,If set to "0", wind speed in the weather file will be read without change.
- `rain_path`: see comments for the `CCfn` function in the `mi.cc.jl` file. If set to "0", rainfall in the weather file will be read without change.
- `site latitude`: latitude of site (decimal degrees)
- `weather station hgt`: height of weather station from ground (m)
- `year`: year of field planting
- `month`: month of field planting (enter 1 for Jan, 2 for Feb, …, 12 for Dec)
- `day`: day of field planting
- `tree age`: age of tree at field planting date (days), usually 365 days (for one year old)
- `planting density`: planting density or no. of palms per hectare (palms/ha)
- `trunk height`: height of the trunk at field planting (m), usually set to 0
- `thinning planting density`: thin to this final planting density (palms/ha) (-1 for no thinning)
- `thinning tree age`: age of tree when thinning was done (days); ignored if `thinning planting density` is set to -1
- `female prob`: probability of getting female flowers (0 - 1), where 0 is never female and 1 is always female, 0.5 is 50-50 chance of female
- `pinnae wgt`: initial pinnae dry weight (kg DM/palm) at time of field planting
- `rachis wgt`: initial rachis dry weight (kg DM/palm) at time of field planting
- `trunk wgt`: initial trunk dry weight (kg DM/palm) at time of field planting
- `roots wgt`: initial roots dry weight (kg DM/palm) at time of field planting
- `maleflo wgt`: initial male flowers dry weight (kg DM/palm) at time of field planting
- `femaflo wgt`: initial female flowers dry weight (kg DM/palm) at time of field planting
- `bunches wgt`: initial bunches dry weight (kg DM/palm) at time of field planting
- `no. of intervals`: no. of subintervals per day for daily soil water flow integration (set to 100 or more time steps in a day)
- `rooting depth`: initial rooting depth (m) at time of field planting
- `max. rate rooting depth`: max. rate of increase in rooting depth (m/day)
- `has watertable`: `true` for water table presence, else `false` for none
- `no. of soil layers`: no. of soil layers (at least 2)
- `layers`: soil layer properties, where each layer comprises four parameters, which are:
    - `thick`: layer thickness (m)
    - `vwc`: initial volumetric soil water content (m3/m3)
    - `clay`: amount of clay in the layer (%)
    - `sand`: amount of sand in the layer (%)
    - note: if you have 3 soil layers, you will need 3 sets of `’thick`, `vwc`, `clay`, and `sand` parameters for each soil layer. See model input file `input.json` for example.
- `sla`: table of tree age vs. specific leaf area (m2 leaf/kg DM leaf); normally use default
- `pinnae n`: table of tree age vs. pinnae N content (fraction); normally use default
- `pinnae m`: table of tree age vs. pinnae mineral content (fraction); normally use default
- `rachis n`: table of tree age vs. rachis N content (fraction); normally use default
- `rachis m`: table of tree age vs. rachis mineral content (fraction); normally use default
- `roots n`: table of tree age vs. roots N content (fraction); normally use default
- `roots m`: table of tree age vs. roots mineral content (fraction); normally use default
- `trunk n`: table of tree age vs. trunk N content (fraction); normally use default
- `trunk m`: table of tree age vs. trunk mineral content (fraction); normally use default

Refer to model input file `input.json` for sample model input entries. You would normally change only a few of these entries, keeping the rest at their default values.

## Daily weather file
Daily weather file for a given site must be prepared. The file must have at least four daily weather parameters: min. and max. air temperature (deg. Celsius), mean daily wind speed (m/s), and daily rainfall amount (mm). Refer to the sample Merlimau weather file `merlimau-wthr.csv`.

Weather data must be arranged by rows to represent the weather for each day, and each column represents the weather properties. All values are separated by comma (comma-delimited values or CSV). The weather file must be in plain text format.

For instance, the first few lines in the sample Merlimau weather file are as follows:

```
year,month,day,tmin,tmax,wind,rain
# Merlimau, Melaka; 2.253213 N, 102.451753 E; 1987/1/1 - 2008/12/31
1987,1,1,23.2,36.6,1.5,0
1987,1,2,23.6,32.9,2.5,0
1987,1,3,21.3,35.1,2,0
1987,1,4,21.9,37.2,2,0
1987,1,5,21.7,36.4,2.1,0
```

where the first row lists the headers so that the `year`, `month`, and `day` are the first three columns for date, followed by the weather parameters `tmin` and `tmax` for min. and max. daily air temperature (deg. C), respectively, `wind` for mean daily wind speed (m/s), and `rain` for amount of daily rainfall (mm/day).

Headers must always be placed on the first line. Please note that these headers are case-sensitive.

The second row (prefixed by the hash character, `#`) is the comment line. The rest of this line will not be read, and its purpose is solely for documentation, such as to describe the weather data, such as site and GPS location and period of weather data. You may have more than one comment line, but all comment lines must be continuous, e.g., from lines 2 to 4 (without a break). Comment lines must always start on the second line.

IMPORTANT: The weather data must be complete, without gaps in data or missing values.

## API
The source code for Sawit.jl is heavily commented. Kindly refer to the code documentation.

## Oil palm data availability
Kindly contact Sime Darby Plantation Research Sdn. Bhd., Malaysia for permission obtaining the historical oil palm data used to validate this model.

## References
Teh, C. B. S., Cheah, S. S., & Kulaveerasingam, H. (2024). Development and validation of an oil palm model for a wide range of planting densities and soil textures in Malaysian growing conditions. Heliyon. [https://doi.org/10.1016/j.heliyon.2024.e32561](https://doi.org/10.1016/j.heliyon.2024.e32561)

Teh, C. B. S. & Cheah, S. S. (2023). Modelling of growth and FFB yield for increasing oil palm productivity. Book of Abstracts. Agriculture, Biotechnology & Sustainability (ABS). MPOB International Palm Oil Congress and Exhibition. Navigating Uncertainties Building Resilience (p. 19). PIPOC 2023. Kajang: Malaysian Palm Oil Board (MPOB) Press. [PIPOC 2023 Book of Abstracts](https://www.researchgate.net/profile/Christopher-Teh-2/publication/376714242_Modelling_of_Growth_and_FFB_Yield_for_Increasing_Oil_Palm_Productivity/links/65845c263c472d2e8e779dc4/Modelling-of-Growth-and-FFB-Yield-for-Increasing-Oil-Palm-Productivity.pdf)

Siang, C. S., Wahid, A. A. A., & Teh, C. B. S. (2022) Standing biomass, dry-matter production and nutrient demand of Tenera oil palm. Agronomy, 12, 426. [https://doi.org/10.3390/agronomy12020426](https://doi.org/10.3390/agronomy12020426)

Cheah, S. S., Teh, C. B. S. (2020). Parameterization of the Farquhar-von Caemmerer-Berry C3 photosynthesis model for oil palm. Photosynthetica, 58(3): 769-779. [https://doi.org/10.32615/ps.2020.020](https://doi.org/10.32615/ps.2020.020)

Cheah, S. S., Teh, C. B. S., Mohd Razi, I., & Mohd Rafii, Y. (2020). Modelling hourly air temperature, relative humidity and solar irradiance over several major oil palm growing areas in Malaysia. Journal of Oil Palm Research, 32(1): 34-49. [https://doi.org/10.21894/jopr.2020.0010](https://doi.org/10.21894/jopr.2020.0010)

Teh, C. B. S., & Cheah, S. S. (2018). Modelling crop growth and yield in palm oil cultivation. In A. Rival (Ed.), Achieving sustainable cultivation of oil palm (Vol. 1, pp. 183-227). Burleigh Dodds Science Publishing. [https://doi.org/10.19103/AS.2017.0018.10](https://doi.org/10.19103/AS.2017.0018.10)

