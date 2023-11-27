"""
    solar_declination(doy::Int)

Solar declination (radians).

# Arguments
- `doy::Int`: day of year

# Returns
- `Float64`: solar declination (radians)
"""
solar_declination(doy::Int) = -0.4093 * cos(0.0172 * (doy + 10))


"""
    solar_constant(doy::Int)

Solar constant (W/m2), corrected for eccentricity.

# Arguments
- `doy::Int`: day of year

# Returns
- `Float64`: solar constant (W/m2)
"""
solar_constant(doy::Int) = 1370 * (1 + 0.033 * cos(0.0172 * (doy - 10)))


"""
    et_solar_radiation(sc, inc)

Instantaneous extra-terrestrial solar irradiance (W/m2).

# Arguments
- `sc`: solar constant (W/m2)
- `inc`: solar inclination (radians)

# Returns
- `Float64`: extra-terrestrial solar irradiance (W/m2)
"""
et_solar_radiation(sc, inc) = max(0, sc * cos(inc))


"""
    saturated_vapor_pressure(ta)

Saturated air vapor pressure (mbar).

# Arguments
- `ta`: air temperature (°C)

# Returns
- `Float64`: saturated air vapor pressure (mbar)
"""
saturated_vapor_pressure(ta) = svp_fn(ta)


"""
    vapor_pressure(ta, tdew)

Instantaneous air vapor pressure (mbar).

# Arguments
- `ta`: air temperature (°C)
- `tdew`: dew temperature (°C)

# Returns
- `Float64`: air vapor pressure (mbar)
"""
vapor_pressure(ta, tdew) = svp_fn(min(ta, tdew))


"""
    vapor_pressure_deficit(svp, vp)

Air vapor pressure deficit (mbar).

# Arguments
- `svp`: saturated air vapor pressure (mbar)
- `vp`: air vapor pressure (mbar)

# Returns
- `Float64`: air vapor pressure deficit (mbar)
"""
vapor_pressure_deficit(svp, vp) = svp - vp


"""
    relative_humidity(svp, vp)

Air relative humidity (RH) (%).

# Arguments
- `svp`: saturated air vapor pressure (mbar)
- `vp`: air vapor pressure (mbar)

# Returns
- `Float64`: air relative humidity (%)
"""
relative_humidity(svp, vp) = 100 * vp / svp


"""
    svp_fn(ta)

Saturated vapor pressure (mbar) at a given air temperature (°C).

# Arguments
- `ta`: air temperature (°C)

# Returns
- `Float64`: saturated vapor pressure (mbar)
"""
svp_fn(ta) = 6.1078 * exp(17.269 * ta / (ta + 237.3))


"""
    _ab(lat, dec)

Internal use.

# Arguments
- `lat`: site latitude (radians)
- `dec`: solar declination (radians)

# Returns
- `NamedTuple{(:a, :b), Tuple{Float64, Float64}}`: `a` and `b` are intermediary variables
"""
function _ab(lat, dec)
    a = sin(lat) * sin(dec)
    b = cos(lat) * cos(dec)
    (; a, b)
end


"""
    twilight_hours(lat, dec)

Local solar time of sunrise and sunset (hours).

# Arguments
- `lat`: site latitude (radians)
- `dec`: solar declination (radians)

# Returns
- `Tuple{Float64, Float64}`: local solar time of sunrise and sunset (hours)
"""
function twilight_hours(lat, dec)
    ab = _ab(lat, dec)
    tss = 12 + (12 / π) * acos(-ab.a / ab.b)  # sunset
    tsr = 24 - tss  # sunrise
    tsr, tss
end


"""
    day_et_solar_radiation(doy::Int, lat, dec)

Extra-terrestrial (outside Earth) solar irradiance (MJ/m2/day).

# Arguments
- `doy::Int`: day of year
- `lat`: site latitude (radians)
- `dec`: solar declination (radians)

# Returns
- `Float64`: extra-terrestrial solar irradiance (MJ/m2/day)
"""
function day_et_solar_radiation(doy::Int, lat, dec)
    ab = _ab(lat, dec)
    aob = ab.a / ab.b
    sc = solar_constant(doy)
    0.027501974 * sc * (ab.a * acos(-aob) + ab.b * sqrt(1 - aob^2))
end


"""
    itg_solar_radiation(lat, mo::Output)

Instantaneous total solar irradiance and its direct and diffuse components (W/m2).

# Arguments
- `lat`: site latitude (radians)
- `mo::Output`: model output

# Notes
- Internal use by function `day_solar_radiation`. This function is called to
  integrate the hourly solar irradiance.

# Returns
- `NamedTuple{Float64, Float64, Float64}(:totrad, :drrad, :dfrad)`:
     total, direct, and diffuse radiation (W/m2).
"""
function itg_solar_radiation(lat, mo::Output)
    @unpack sunrise, sunset, doy, tmin, tmax, solardec = mo.op
    solarcon = solar_constant(doy)
    dew_temp = 23.0

    function hour_wthr(th)  # th = local solar time (hours)
        air_temp = air_temperature(th, tmin, tmax, sunrise, sunset)
        svp = saturated_vapor_pressure(air_temp)
        vp = vapor_pressure(air_temp, dew_temp)
        rh = relative_humidity(svp, vp)
        solarinc, _, _ = solar_position(th, lat, solardec)
        etrad = et_solar_radiation(solarcon, solarinc)
        totrad, drrad, dfrad = solar_radiation(rh, etrad)
        (; totrad, drrad, dfrad)
    end
end


"""
    day_solar_radiation(lat, mo::Output)

Daily total solar irradiance and its direct and diffuse components (MJ/m2/day).

# Arguments
- `lat`: site latitude (radians)
- `mo::Output`: model output

Returns:
- `Tuple{Float64, Float64, Float64}`:
        daily total, direct, and diffuse solar irradiance (MJ/m2/day)
"""
function day_solar_radiation(lat, mo::Output)
    op = mo.op
    tsr, tss = op.sunrise, op.sunset
    fields = [:totrad, :drrad, :dfrad]   # solar radiation components
    hour_rad = itg_solar_radiation(lat, mo)
    itg = integrate(5, tsr, tss)
    res = itg(hour_rad; fields=fields)
    res .* 3600 / 10^6  # in MJ/m2/day
end


"""
    solar_position(th, lat, dec)

Solar position: inclination, height (solar elevation), and azimuth (radians).

# Arguments
- `th`: local solar time (hours)
- `lat`: site latitude (radians)
- `dec`: solar declination (radians)

# Notes
- Valid only from sunrise to sunset hours.

# Returns
- `Tuple{Float64, Float64, Float64}`: solar inclination, height, and azimuth (radians)
"""
function solar_position(th, lat, dec)
    ab = _ab(lat, dec)
    ha = π / 12 * (th - 12)  # hour angle

    inc = min(0.5 * π, acos(ab.a + ab.b * cos(ha)))  # solar inclination
    hgt = 0.5 * π - inc  # solar height/elevation

    # azimuth (angle from North in a clockwise direction):
    a = cos(dec) * (cos(lat) * tan(dec) + sin(lat) * cos(ha)) / cos(hgt)
    acosa = acos(max(-1, min(1, a)))
    azi = th <= 12 ? acosa : (2π - acosa)

    inc, hgt, azi
end


"""
    solar_radiation(rh, etrad)

Instantaneous total solar irradiance and its direct and diffuse components (W/m2).

# Arguments
- `rh`: relative humidity (%)
- `etrad`: ET solar irradiance (W/m2)

# Returns
- `Tuple{Float64, Float64, Float64}`:
        instantaneous total, direct, and diffuse solar irradiance (W/m2)
"""
function solar_radiation(rh, etrad)
    kt = 1.1595 - 0.0106 * rh  # sky clearness index = total radiation / ET radiation
    it = etrad * kt    # total
    # diffuse partitioning based on Khatib et al. (2012) for 28 Malaysian sites:
    ratio = 0.9505 + 0.91634 * kt - 4.851 * kt^2 + 3.2353 * kt^3
    idf = it * ratio    # diffuse
    idr = it - idf      # direct
    it, idr, idf
end


"""
    net_solar_radiation(ta, vp, totrad)

Instantaneous net solar irradiance (W/m2).

# Arguments
- `ta`: air temperature (°C)
- `vp`: air vapor pressure (mbar)
- `totrad`: total solar irradiance (W/m2)

# Returns
- `Float64`: instantaneous net solar irradiance (W/m2)
"""
function net_solar_radiation(ta, vp, totrad)
    p = 0.15  # reflection
    sb = 5.67 * 10^-8  # Stefan-Boltzmann constant (W/m2/K4)
    tak = ta + 273.15  # air temperature in Kelvin (K)
    rnl = 0.98 * sb * tak^4 * (1.31 * (vp / tak)^(1 / 7) - 1)  # net longwave
    (1 - p) * totrad + rnl
end


"""
    air_temperature(th, tmin, tmax, tsr, tss, lag=1.5)

Instantaneous air temperature (°C).

# Arguments
- `th`: local solar time (hours)
- `tmin`: min. air temperature (°C)
- `tmax`: max. air temperature (°C)
- `tsr`: local solar time of sunrise (hours)
- `tss`: local solar time of sunset (hours)

# Keywords
- `lag`: no. of hours after sunrise that air temperature starts to increase (hours)

# Returns
- `Float64`: instantaneous air temperature (°C)
"""
function air_temperature(th, tmin, tmax, tsr, tss, lag=1.5)
    dl = tss - tsr   # day length
    if (tsr + lag) <= th <= tss
        n1 = π * (th - tsr - lag) / dl
        ta = tmin + (tmax - tmin) * sin(n1)
    else
        tset = tmin + (tmax - tmin) * sin(π * (dl - lag) / dl)
        t = th < (tsr + lag) ? tsr : -tss
        ta = tset + ((tmin - tset) * (th + t)) / ((tsr + lag) + tsr)
    end
    ta
end


"""
    svp_slope(ta)

Slope of the curve between saturated air vapor pressure (SVP) and
air temperature (mbar/°C).

# Arguments
- `ta`: air temperature (°C)

# Returns
- `Float64`: slope of SVP vs. air temperature curve (mbar/°C)
"""
function svp_slope_fn(ta)
    n1 = exp(17.269 * ta / (ta + 237.3))
    n2 = (ta + 237.3) ^ 2
    25029.4 * n1 / n2
end


"""
    wind_speed(th, uday, tsr, tss, lag=1.5)

Instantaneous wind speed (m/s).

# Arguments
- `th`: local solar time (hours)
- `uday`: daily mean wind speed (m/s)
- `tsr`: local solar time of sunrise (hours)
- `tss`: local solar time of sunset (hours)

# Keywords
- `lag`: no. of hours after sunrise that wind speed starts to increase (hours)

# Returns
- `Float64`: instantaneous wind speed (m/s)
"""
function wind_speed(th, uday, tsr, tss, lag=1.5)
    umin = 0.559134814 * uday^1.25  # min. wind speed
    umax = 1.797613613 * uday^0.75  # max. wind speed
    # wind speed varies only between (tsr+lag) and (tss+lag) hours:
    dl = tss - tsr    # day length
    Δ = (umax - umin) * sin(π * (th - tsr - lag) / dl)
    bwithin = (tsr + lag) <= th <= (tss + lag)
    bwithin ? (umin + Δ) : umin
end


"""
    read_wthr(date::Date, wthr::AbstractDataFrame)

Read in daily weather data for the given date.

# Arguments
- `date::Date`: date of weather retrieval
- `wthr::AbstractDataFrame`: full weather data

# Returns
- `Tuple{Float64, Float64, Float64, Float64}`:
        daily min. and max. air temperature, wind speed, and rain
"""
function read_wthr(date::Date, wthr::AbstractDataFrame)
    df0 = @view wthr[wthr.date .== date, :]
    tmin = df0.tmin[1]
    tmax = df0.tmax[1]
    wind = df0.wind[1]
    rain = df0.rain[1]
    tmin, tmax, wind, rain
end


"""
    est_hour_wthr(lat, mo::Output)

Estimate the hourly meteorological paramaters based on given daily parameters.

# Arguments
- `lat`: site latitude (radians)
- `mo::Output`: model output

# Notes
- Does not determine hourly rainfall.

# Returns
- `Function`: Returns a function object that accepts a single argument for local
              solar time (hours) to estimate the hourly weather/meterological
              properties (except rain). This function object returns `OPHour`.
"""
function est_hour_wthr(lat, mo::Output)
    @unpack sunrise, sunset, doy, date, tmin, tmax, solardec, wind = mo.op
    dew_temp = 23.0
    solarcon = solar_constant(doy)

    function hour_wthr(th)  # th = local solar time (hours)
        air_temp = air_temperature(th, tmin, tmax, sunrise, sunset)
        svp = saturated_vapor_pressure(air_temp)
        vp = vapor_pressure(air_temp, dew_temp)
        rh = relative_humidity(svp, vp)
        vpd = vapor_pressure_deficit(svp, vp)
        svp_slope = svp_slope_fn(air_temp)
        solarinc, solarhgt, solarazi = solar_position(th, lat, solardec)
        etrad = et_solar_radiation(solarcon, solarinc)
        totrad, drrad, dfrad = solar_radiation(rh, etrad)
        netrad = net_solar_radiation(air_temp, vp, totrad)
        hrwind = wind_speed(th, wind, sunrise, sunset)  # `wind` = daily wind speed

        OPHour(
            date = date,
            doy = doy,
            solarhour = th,
            air_temp = air_temp,
            dew_temp = dew_temp,
            svp = svp,
            vp = vp,
            rh = rh,
            vpd = vpd,
            svp_slope = svp_slope,
            solarinc = solarinc,
            solarhgt = solarhgt,
            solarazi = solarazi,
            solarcon = solarcon,
            etrad = etrad,
            totrad = totrad,
            drrad = drrad,
            dfrad = dfrad,
            netrad = netrad,
            wind = hrwind
        )
    end
end


"""
    update_wthr!(mo::Output, mi::Input)

Update the daily meteorological paramaters.

# Arguments
- `mo::Output`: model output
- `mi::Input`: model input

# Returns
- nothing
"""
function update_wthr!(mo::Output, mi::Input)
    op = mo.op
    date = op.date
    lat = mi.lat
    cc = mi.ccfn

    tmin, tmax, wind, rain = read_wthr(date, mi.wthr)
    op.tmin = cc.Δta(date, tmin)
    op.tmax = cc.Δta(date, tmax)
    op.wind = cc.Δwind(date, wind)
    op.rain = cc.Δrain(date, rain)

    op.doy = doy = dayofyear(date)
    op.solardec = solardec = solar_declination(doy)
    op.sunrise, op.sunset = twilight_hours(lat, solardec)
    op.etrad = day_et_solar_radiation(doy, lat, solardec)
    op.totrad, op.drrad, op.dfrad = day_solar_radiation(lat, mo)
end
