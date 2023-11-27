"""
    CCfn{C<:Function, T<:Function, W<:Function, R<:Function}

Collection of functions to change daily meteorological properties due to climate change.

# Fields
- `co2::C<:Function`: function to determine the ambient CO2 level
- `co2_path::String`: code to denote which CO2 function to use (default: "PAST")
- `Δta::T<:Function`: function to determine the change in air temperature
- `Δta_path::String`: code to denote which temperature function to use (default: "0")
- `Δwind::W<:Function`: function to determine the change in wind speed
- `Δwind_path::String`: code to denote which wind speed function to use (default: "0")
- `Δrain::R<:Function`: function to determine the change in rainfall amount
- `Δrain_path::String`: code to denote which rain function to use (default: "0")
"""
@with_kw struct CCfn{C<:Function, T<:Function, W<:Function, R<:Function}
    co2::C = co2_past
    co2_path::String = "PAST"

    Δta::T = Δta_const(0.0)
    Δta_path::String = "0"

    Δwind::W = Δwind_const(0.0)
    Δwind_path::String = "0"

    Δrain::R = Δrain_const(0.0)
    Δrain_path::String = "0"
end


"""
     CCfn(dct::AbstractDict)

Construct CCfn object from the model input.

# Arguments
- `dct::AbstractDict`: model input stored in a dictonary

# Notes
- model input dictionary `dct` must have the following keys as `String`:

     "co2_path", "ta_path", "wind_path", and "rain_path"

  where `co2_path` can be set to (case insensitive):
      - blank or "past" to estimate CO2 levels from past CO2 trends
      - "rcp2.6" for RCP2.6 path
      - "rcp4.5" for RCP4.5 path
      - "rcp6" for RCP6 path
      - "rcp8.5" for RCP8.5 path
      - a constant absolute level, e.g., "450" or "550" for 450 and 550
        µmol CO2/mol air, respectively.

  where `ta_path` can be set to (case insensitive):
      - blank for no change (+0 °C)
      - "rcp2.6" for RCP2.6 path
      - "rcp4.5" for RCP4.5 path
      - "rcp6" for RCP6 path
      - "rcp8.5" for RCP8.5 path
      - a constant change, e.g., "1" or "2" for an increase by +1 and +2 °C, respectively

  where `wind_path` can be set to (case insensitive):
      - blank for no change (+0 m/s)
      - a constant absolute change, e.g., "1" or "-0.5" for an increase by
        +1 m/s and a decrease by -0.5 m/s, respectively.

  where `rain_path` can be set to (case insensitive):
      - blank for no change (+0%)
      - a constant percent change, e.g., "10" or "-20" for an increase by
        10% and a decrease by -20%, respectively.

  where all constants, even though they are numbers, must be still be enclosed by
  double quotes, as they will be read as `String`.

# Returns
- `CCFn`: new object instantiated
"""
function CCfn(dct::AbstractDict)
    cleanup(txt) = uppercase(strip(txt))
    to_val(txt) = isempty(txt) ? 0.0 : parse(Float64, txt)

    co2_path = cleanup(dct["co2_path"])
    if isempty(co2_path) || co2_path == "PAST"
        co2 = co2_past
    elseif co2_path == "RCP2.6"
        co2 = co2_rcp2_6
    elseif co2_path == "RCP4.5"
        co2 = co2_rcp4_5
    elseif co2_path == "RCP6"
        co2 = co2_rcp6
    elseif co2_path == "RCP8.5"
        co2 = co2_rcp8_5
    else
        co2 = co2_const(to_val(co2_path))
    end

    Δta_path = cleanup(dct["ta_path"])
    if Δta_path == "RCP2.6"
        Δta = Δta_rcp2_6
    elseif Δta_path == "RCP4.5"
        Δta = Δta_rcp4_5
    elseif Δta_path == "RCP6"
        Δta = Δta_rcp6
    elseif Δta_path == "RCP8.5"
        Δta = Δta_rcp8_5
    else
        Δta = Δta_const(to_val(Δta_path))
    end

    Δwind_path = cleanup(dct["wind_path"])
    Δwind = Δwind_const(to_val(Δwind_path))

    Δrain_path = cleanup(dct["rain_path"])
    Δrain = Δrain_const(to_val(Δrain_path))

    CCfn(
        co2 = co2,
        co2_path = co2_path,
        Δta = Δta,
        Δta_path = Δta_path,
        Δwind = Δwind,
        Δwind_path = Δwind_path,
        Δrain = Δrain,
        Δrain_path = Δrain_path,
    )
end


"""
    co2_past(date::Date)

Ambient CO2 level (µmol CO2/mol air), estimated from actual/measured past trends.

# Arguments
- `date::Date`: date at which to estimate the CO2 level

# Returns
- `Float64`: ambient CO2 level (µmol CO2/mol air)
"""
function co2_past(date::Date)
    y = fraction_year(date)
    sqrt(39413600 - 40620.1096 * y + 10.49094 * y^2)
end


"""
    co2_rcp2_6(date::Date)

Ambient CO2 level (µmol CO2/mol air) based on RCP2.6 pathway (best case scenario).

# Arguments
- `date::Date`: date at which to estimate the CO2 level

# Returns
- `Float64`: ambient CO2 level (µmol CO2/mol air)
"""
function co2_rcp2_6(date::Date)
    y = fraction_year(date)
    a, b, c, d, e = 18.98954528, -0.00099049, -0.01883991, 2.4541E-07, 4.67579E-06
    ((a + c * y + e * y^2) / (1 + b * y + d * y^2)) ^ 2
end


"""
    co2_rcp4_5(date::Date)

Ambient CO2 level (µmol CO2/mol air) based on RCP4.5 pathway (most likely scenario).

# Arguments
- `date::Date`: date at which to estimate the CO2 level

# Returns
- `Float64`: ambient CO2 level (µmol CO2/mol air)
"""
function co2_rcp4_5(date::Date)
    y = fraction_year(date)
    a, b, c, d, e = 440.1388926, -0.0443354, -19.5816588, 0.000491542, 0.217853519
    (a + c * sqrt(y) + e * y) / (1 + b * sqrt(y) + d * y)
end


"""
    co2_rcp6(date::Date)

Ambient CO2 level (µmol CO2/mol air) based on RCP6 pathway.

# Arguments
- `date::Date`: date at which to estimate the CO2 level

# Returns
- `Float64`: ambient CO2 level (µmol CO2/mol air)
"""
function co2_rcp6(date::Date)
    y = fraction_year(date)
    (-0.01062159 + 53339.74186 / y^2) ^ -1
end


"""
    co2_rcp8_5(date::Date)

Ambient CO2 level (µmol CO2/mol air) based on RCP8.5 pathway (BAU scenario).

# Arguments
- `date::Date`: date at which to estimate the CO2 level

# Returns
- `Float64`: ambient CO2 level (µmol CO2/mol air)
"""
function co2_rcp8_5(date::Date)
    y = fraction_year(date)
    (-0.01513419 + 71336.45572 / y^2) ^ -1
end


"""
    co2_const(val)

Constant/fixed ambient CO2 level (µmol CO2/mol air).

# Arguments
- `val`: desired CO2 level (µmol CO2/mol air)

# Notes
- This function returns an anonymous function. Unlike its counterparts,
  e.g., `co2_past` and `co2_rcp4_5`, that accepts a date to estimate
  the CO2 level, this function ignores the passed argument of date and
  returns only a pre-defined constant value.

# Returns
- `Float64`: ambient CO2 level (µmol CO2/mol air)
"""
co2_const(val) = (_ -> val)


"""
    Δta_rcp2_6(date::Date)

Ambient air temperature, after increase following the RCP2.6 pathway (best case).
Based on 3 °C climate sensitivity.

# Arguments
- `date::Date`: date at which to estimate the air temperature
- `ta`: daily air temperature (°C)

# Returns
- `Float64`: daily air temperature after increase (°C)
"""
function Δta_rcp2_6(date::Date, ta)
    y = fraction_year(date)
    a, b, c, d = -6275.67818, 9.045301705, -0.00434475, 6.95626E-07
    ta + (a + b * y + c * y^2 + d * y^3) ^ 2
end


"""
    Δta_rcp4_5(date::Date)

Ambient air temperature, after increase following the RCP4.5 pathway (most likely).
Based on 3 °C climate sensitivity.

# Arguments
- `date::Date`: date at which to estimate the air temperature
- `ta`: daily air temperature (°C)

# Returns
- `Float64`: daily air temperature after increase (°C)
"""
function Δta_rcp4_5(date::Date, ta)
    y = fraction_year(date)
    a, b, c = 8.963069745, -0.00057523, -0.00454612
    ta + (a + c * y) / (1 + b * y)
end


"""
    Δta_rcp6(date::Date)

Ambient air temperature, after increase following the RCP6 pathway.
Based on 3 °C climate sensitivity.

# Arguments
- `date::Date`: date at which to estimate the air temperature
- `ta`: daily air temperature (°C)

# Returns
- `Float64`: daily air temperature after increase (°C)
"""
function Δta_rcp6(date::Date, ta)
    y = fraction_year(date)
    a, b, c = 239.6537001, -0.22571417, 5.32018E-05
    ta + (a + b * y + c * y^2) ^ -1
end


"""
    Δta_rcp8_5(date::Date)

Ambient air temperature, after increase following the RCP8.5 pathway (BAU).
Based on 3 °C climate sensitivity.

# Arguments
- `date::Date`: date at which to estimate the air temperature
- `ta`: daily air temperature (°C)

# Returns
- `Float64`: daily air temperature after increase (°C)
"""
function Δta_rcp8_5(date::Date, ta)
    y = fraction_year(date)
    a, b, c = 6868.386721, -6.88774894, 0.00172699
    ta + sqrt(a + b * y + c * y^2)
end


"""
    Δta_const(val)

Change ambient air temperature by a constant amount.

# Arguments
- `val`: amount in change (°C); negative value for decrease

# Notes
- This function returns an anonymous function that adds a
  constant value to the air temperature (`ta`).

# Returns
- `Float64`: modified air temperature (°C)
"""
Δta_const(val) = ((_, ta) -> ta + val)


"""
    Δwind_const(val)

Change wind speed by a constant amount.

# Arguments
- `val`: amount in change (m/s); negative value for decrease

# Notes
- This function returns an anonymous function that adds a
  constant value to the wind speed (`wind`).

# Returns
- `Float64`: modified wind speed (m/s)
"""
Δwind_const(val) = ((_, wind) -> wind + val)


"""
    Δrain_const(pct)

Change rainfall amount by a constant percentage (not absolute value).

# Arguments
- `pct`: percentage change (%); negative % for decrease

# Notes
- This function returns an anonymous function that changes the
  rainfall (`rain`) by a constant percentage.

# Returns
- `Float64`: modified rainfall amount (mm)
"""
Δrain_const(pct) = ((_, rain) -> rain * (1 + 0.01 * pct))
