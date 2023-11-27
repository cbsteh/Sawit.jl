module GoF

using Statistics


# NOTE:
# Function arguments always take in observed, then estimated values.
# Differences between model estimations and measured/observed values
# are always calculated as: estimations - observations.

"""
    filter_data(obs::AbstractVector, est::AbstractVector)

Remove zero and missing values from observations and their corresponding
estimated values.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Some GoF indexes will have divide-by-zero error if there are 0 values
  in observations. Use this function to remove all 0 from observations,
  and their corresponding estimated values.

# Returns
- `Tuple{AbstractVector, AbstractVector}`: observations and estimations
"""
function filter_data(obs::AbstractVector, est::AbstractVector)
    idxs = findall(v -> !ismissing(v) && !isapprox(v, 0), obs)
    o = getindex(obs, idxs)
    e = getindex(est, idxs)
    (; o, e)
end


"""
    minmax(obs::AbstractVector, est::AbstractVector)

Minimum and maximum difference (errors) between estimated and observed values.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `NamedTuple{(:mn, :mx), {Float64, Float64}}`: smallest and largest errors
"""
function minmax(obs::AbstractVector, est::AbstractVector)
    mn, mx = extrema(est .- obs)
    (; mn, mx)
end


"""
    quartiles(obs::AbstractVector, est::AbstractVector)

First, second, and third quartiles of errors (difference between estimated
and observed values).

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `NamedTuple{(:Q1, :Q2, :Q3), {Float64, Float64, Float64}}`:
     quartile 1 (25th-percentile), 2 (median), and 3 (75th-percentile) of errors
"""
function quartiles(obs::AbstractVector, est::AbstractVector)
    err = est .- obs
    (Q1=quantile(err, 0.25), Q2=quantile(err, 0.5), Q3=quantile(err, 0.75))
end


"""
    NMAE(obs::AbstractVector, est::AbstractVector)

Normalized mean absolute error.

0 <= NMAE < +INF. Best fit = 0. Large +ve = large overall error.
Example: NMAE = 0.02 means an average absolute difference by 2% of the observations.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: normalized mean absolute error
"""
function NMAE(obs::AbstractVector, est::AbstractVector)
    (est .- obs .|> abs |> sum) / sum(obs)
end


"""
    MAE(obs::AbstractVector, est::AbstractVector)

Mean absolute error.

0 <= MAE < +INF. Best fit = 0.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: mean absolute error
"""
function MAE(obs::AbstractVector, est::AbstractVector)
    est .- obs .|> abs |> mean
end


"""
    NMBE(obs::AbstractVector, est::AbstractVector)

Normalized mean bias error.

-INF < NMBE < +INF. Best fit = 0, which means no overall bias.
Large +ve = overestimate, large -ve = underestimate.
Example: NMBE = +0.02 means an average overestimation by 2% of the observations.

# Notes
- Errors can cancel out, i.e., positve errors cancel out negative
  to give small or even zero net error.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: normalized mean bias error
"""
function NMBE(obs::AbstractVector, est::AbstractVector)
    sum(est .- obs) / sum(obs)
end


"""
    MBE(obs::AbstractVector, est::AbstractVector)

Mean bias error.

-INF < MBE < +INF. Best fit = 0.
Large +ve = overestimate, large -ve = underestimate.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: mean bias error
"""
function MBE(obs::AbstractVector, est::AbstractVector)
    mean(est .- obs)
end


"""
    VE(obs::AbstractVector, est::AbstractVector)

Volumetric efficiency.

0<= VE <= 1. Best fit = 1

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Criss, R. E., & Winston, W. E. (2008), Do Nash values have value? Discussion and
  alternate proposals. Hydrological Processes. 22, 2723-2725.

# Returns
- `Float64`: percent bias (%)
"""
function VE(obs::AbstractVector, est::AbstractVector)
    1 - (est .- obs .|> abs |> sum) / sum(obs)
end


"""
    PBIAS(obs::AbstractVector, est::AbstractVector)

Percent bias (%): mean tendency of the model to over- or underestimate.

-INF <= PBIAS < +INF. Best fit = 0%.
Positive values indicate overestimation, whereas negative underestimation.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Yapo P. O., Gupta H. V., & Sorooshian S. (1996). Automatic calibration of
  conceptual rainfall-runoff models: sensitivity to calibration data.
  Journal of Hydrology, 181(i1-4), 23-48.
- Sorooshian, S., Duan, Q., & Gupta. V. K. (1993). Calibration of rainfall-runoff
  models: Application of global optimization to the Sacramento Soil Moisture
  Accounting Model. Water Resources Research, 29(4), 1185-1194.

# Returns
- `Float64`: percent bias (%)
"""
function PBIAS(obs::AbstractVector, est::AbstractVector)
    100 * (sum(est .- obs) / sum(obs))
end


"""
    MAPE(obs::AbstractVector, est::AbstractVector)

Mean absolute percentage error (%).

0% <= MAPE < +INF. Best fit = 0%.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Fails if there are 0 values in observations.

# Returns
- `Float64`: mean absolute percentage error (%)
"""
function MAPE(obs::AbstractVector, est::AbstractVector)
    100 * ((est .- obs) ./ obs .|> abs |> mean)
end


"""
    MBPE(obs::AbstractVector, est::AbstractVector)

Mean bias percentage error (%).

-INF < MBPE < +INF. Best fit = 0%.
Large +ve = overestimate, large -ve = underestimate.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Fails if there are 0 values in observations.

# Returns
- `Float64`: mean bias percentage error (%)
"""
function MBPE(obs::AbstractVector, est::AbstractVector)
    100 * ((est .- obs) ./ obs |> mean)
end


"""
    MAAPE(obs::AbstractVector, est::AbstractVector)

Mean arctangent (inverse tangent) absolute percentage error (radians).

0 <= MAAPE <= π/2. Best fit = 0 radians, Worst fit = π/2 (i.e.,
max. angle or observed are perpendicular to estimated values).

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Fails if there are 0 values in observations.
- Kim, S., & Kim, H. (2016). A new metric of absolute percentage error for
  intermittent demand forecasts. International Journal of Forecasting, 32(3), 669-679.

# Returns
- `Float64`: mean arctangent absolute percentage error (radians)
"""
function MAAPE(obs::AbstractVector, est::AbstractVector)
    (est .- obs) ./ obs .|> abs .|> atan |> mean
end


"""
    MDAPE(obs::AbstractVector, est::AbstractVector)

Median absolute percentage error (%).

0% <= MDAPE < +INF. Best fit = 0%.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Fails if there are 0 values in observations.
- Like MAPE except median is used instead of mean to reduce any outliers' influence.

# Returns
- `Float64`: median absolute percentage error (%)
"""
function MDAPE(obs::AbstractVector, est::AbstractVector)
    100 * ((est .- obs) ./ obs .|> abs |> median)
end


"""
    RMSE(obs::AbstractVector, est::AbstractVector)

Root mean square error.

0 <= RMSE < +INF. Best fit = 0.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: root mean square error
"""
function RMSE(obs::AbstractVector, est::AbstractVector)
    (est .- obs) .^ 2 |> mean |> sqrt
end


"""
    RSR(obs::AbstractVector, est::AbstractVector)

Ratio of root mean square error (RMSE) to standard deviation.

0 <= RSR < +INF. Best fit = 0.
RSR < 0.5 for very good; 0.5-0.6 for good, 0.6-0.7 for satisfactory.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: ratio of RMSE to standard deviation
"""
function RSR(obs::AbstractVector, est::AbstractVector)
    μo = mean(obs)
    n1 = sum((est .- obs) .^ 2)
    n2 = sum((obs .- μo) .^ 2)
    sqrt(n1) / sqrt(n2)
end


"""
    NSE(obs::AbstractVector, est::AbstractVector)

Nash-Sutcliffe Efficiency index.

-INF < NSE <= 1. Best fit = 1.
NSE > 0.75 for very good; 0.75-0.65 for good and 0.65-0.50 for satisfactory.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Nash, J. E., & Sutcliffe, J. V. (1970). River flow forecasting through
  conceptual models part I — A discussion of principles. Journal of Hydrology,
  10, 282–290.

# Returns
- `Float64`: Nash-Sutcliffe Efficiency index
"""
function NSE(obs::AbstractVector, est::AbstractVector)
    μo = mean(obs)
    n1 = sum((est .- obs) .^ 2)
    n2 = sum((obs .- μo) .^ 2)
    1 - n1 / n2
end


"""
    RNSE(obs::AbstractVector, est::AbstractVector)

Relative Nash-Sutcliffe Efficiency index.

0 <= RNSE <= 1. Best fit = 1.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Nash, J. E., & Sutcliffe, J. V. (1970). River flow forecasting through
  conceptual models part I — A discussion of principles. Journal of Hydrology,
  10, 282–290.

# Returns
- `Float64`: Nash-Sutcliffe Efficiency index
"""
function RNSE(obs::AbstractVector, est::AbstractVector)
    μo = mean(obs)
    n1 = sum(((est .- obs) ./ obs) .^ 2)
    n2 = sum(((obs .- μo) ./ μo) .^ 2)
    1 - n1 / n2
end


"""
    MNSE(obs::AbstractVector, est::AbstractVector)

Modified Nash-Sutcliffe Efficiency index.

0 <= MNSE <= 1. Best fit = 1.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- This original NSE is modified, so that the errors are not inflated by the
  squared values of differences, because the squares are replaced by absolute values.
- Krause, P., Boyle, D. P., & Base, F. (2205). Comparison of different efficiency
  criteria for hydrological model assessment, Advances in Geosciences. 5, 89-97.
- Legates, D. R., & McCabe Jr., G. J. (1999), Evaluating the use of
  "Goodness-of-Fit" measures in hydrologic and hydroclimatic model validation.
  Water Resources Research. 35(1), 233-241.
- Nash, J. E., & Sutcliffe, J. V. (1970). River flow forecasting through
  conceptual models part I — A discussion of principles. Journal of Hydrology,
  10, 282–290.

# Returns
- `Float64`: Nash-Sutcliffe Efficiency index
"""
function MNSE(obs::AbstractVector, est::AbstractVector)
    μo = mean(obs)
    n1 = est .- obs .|> abs |> sum
    n2 = obs .- μo .|> abs |> sum
    1 - n1 / n2
end


"""
    NMSE(obs::AbstractVector, est::AbstractVector)

Normalized mean square error.

-INF < NMSE < +INF. Best fit = 0.
NMSE between -0.5 and +0.5 is considered acceptable.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: normalized mean square error
"""
function NMSE(obs::AbstractVector, est::AbstractVector)
    mean((est .- obs) .^ 2) / (mean(obs) * mean(est))
end


"""
    FB(obs::AbstractVector, est::AbstractVector)

Fractional bias.

-2 <= FB <= +2. Best fit = 0.
FB between -0.5 and +0.5 is considered acceptable.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: fractional bias
"""
function FB(obs::AbstractVector, est::AbstractVector)
    2 * mean((est .- obs) ./ (est .+ obs))
end


"""
    COE(obs::AbstractVector, est::AbstractVector)

Coefficient of Efficiency index.

-INF < COE <= 1. Best fit = 1.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Returns
- `Float64`: fractional bias
"""
function COE(obs::AbstractVector, est::AbstractVector)
    μo = mean(obs)
    n1 = est .- obs .|> abs |> sum
    n2 = obs .- μo .|> abs |> sum
    1 - n1 / n2
end


"""
    MIELKE(obs::AbstractVector, est::AbstractVector)

Revised Mielke index (reduction in r due to model errors).

-1 <= MIELKE <= 1. Best fit = 1.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Duveiller, G., Fasbender, D., & Meroni, M. (2016). Revisiting the concept of
  a symmetric index of agreement for continuous datasets. Scientific Reports,
  6(19401), 1-14.
- Mielke, P. (1984). Meteorological applications of permutation techniques
  based on distance functions. In Krishnaiah, P. & Sen, P. (eds.). Handbook of
  Statistics Vol. 4, 813–830, Elsevier, Amsterdam, The Netherlands.

# Returns
- `Float64`: revised Mielke index
"""
function MIELKE(obs::AbstractVector, est::AbstractVector)
    sdo = std(obs)
    sde = std(est)
    μo = mean(obs)
    μe = mean(est)
    n1 = sdo / sde
    n2 = sde / sdo
    alpha = n1 + n2  + ((μo - μe) .^ 2) / (sdo * sde)
    2 / alpha * cor(obs, est)
end


"""
    PI(obs::AbstractVector, est::AbstractVector)

Persistence index.

-INF < PI <= 1. Best fit = 1, >0 satisfactory, <=0 poor.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Kitanidis, P. K., & Bras, R. L. (1980) Real-time forecasting with a
  conceptual hydrologic model. 2: application and results. Water Resources
  Research, 16(6), 1034–1044.

# Returns
- `Float64`: persistence index
"""
function PI(obs::AbstractVector, est::AbstractVector)
    n1, n2 = 0, 0
    sz = length(obs) - 1
    for i ∈ 1:sz
        n1 += (est[i+1] - obs[i+1])^2
        n2 += (obs[i+1] - obs[i])^2
    end
    1 - n1 / n2
end


"""
    THEIL_U2(naive::AbstractVector, est::AbstractVector)

Theil's coefficient of inequality (U2, 2nd version).

0 <= U2 < +INF. U2 < 1 = is better, 1 = is no better/worse,
>1 = is worse than naive (guess) estimates.

# Arguments
- `naive::AbstractVector`: vector of naive (guess) values
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Thiel, H. (1966). Applied Economic Forecasting. Chicago, Rand McNally.
- Theil, H. (1958). Economic Forecasts and Policy. Amsterdam, North Holland.

# Returns
- `Float64`: Theil's coefficient of inequality (2nd version)
"""
function THEIL_U2(naive::AbstractVector, est::AbstractVector)
    n1, n2 = 0, 0
    sz = length(naive) - 1
    for i ∈ 1:sz
        n1 += ((est[i+1] - naive[i+1]) / naive[i])^2
        n2 += ((naive[i+1] - naive[i]) / naive[i])^2
    end
    sqrt(n1 / n2)
end


"""
    KGE(obs::AbstractVector, est::AbstractVector)

Revised Kling-Gupta Efficiency (KGE) index.

-INF < KGE <= +1. Best fit = 1. KGE should be at least > -0.41
(note: KGE = -0.41 means estimations are all constant and equal to
the observed mean, so that estimations have zero variation and zero
correlation with measurements).

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Based on Kling et al. (2012) which uses CV instead of SD for variability ratio to avoid
  cross correlation between the bias and variability ratio.
- Kling, H., Fuchs, M., & Paulin, M. (2012). Runoff conditions in the upper Danube
  basin under an ensemble of climate change scenarios. Journal of Hydrology, 424, 264-277.
- Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). Decomposition of the
  mean squared error and NSE performance criteria: Implications for improving hydrological
  modelling. Journal of Hydrology, 377, 80–91.

# Returns
- `Float64`: Kling-Gupta Efficiency index
"""
function KGE(obs::AbstractVector, est::AbstractVector)
    r = cor(obs, est)
    μo = mean(obs)
    μe= mean(est)
    cv_obs = std(obs) / μo    # use CV, instead of SD
    cv_est = std(est) / μe
    alpha = (r - 1)^2 + (μe / μo - 1)^2 + (cv_est / cv_obs - 1)^2
    1 - sqrt(alpha)
end


"""
    d(obs::AbstractVector, est::AbstractVector)

Willmott's original index of agreement.

0 <= d <= +1. Best fit = 1.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Willmott, C. J. (1981). On the validation of models. Physical Geography, 2, 184–194.

# Returns
- `Float64`: Willmott's revised index of agreement
"""
function d(obs::AbstractVector, est::AbstractVector)
    μo = mean(obs)
    n1 = est .- obs .|> abs |> sum
    n2 = abs.(est .- μo) .+ abs.(obs .- μo) |> sum
    1 - n1 / n2
end


"""
    dr(obs::AbstractVector, est::AbstractVector)

Willmott's revised index of agreement.

-1 <= dr <= +1. Best fit = 1. Worst fit = -1 (perhaps due to lack of data/variation)

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Willmott, C. J., Robeson, S. M., & Matsuura, K. (2012). A refined index of model
  performance, International Journal of Climatolology, 32, 2088-2094.
- Willmott, C. J. (1981). On the validation of models. Physical Geography, 2, 184–194.

# Returns
- `Float64`: Willmott's revised index of agreement
"""
function dr(obs::AbstractVector, est::AbstractVector)
    μo = mean(obs)
    n1 = est .- obs .|> abs |> sum
    n2 = 2 * (obs .- μo .|> abs |> sum)
    (n1 <= n2) ? (1 - n1 / n2) : (n2 / n1 - 1)
end


"""
    LINCCC(obs::AbstractVector, est::AbstractVector)

Lin's concordance correlation coefficient.

-1 <= LINCCC <= +1. LINCCC = 1 means strong concordance (agreement),
-1 strong discordance, and 0 no concordance.

# Arguments
- `obs::AbstractVector`: vector of observations
- `est::AbstractVector`: vector of estimations (simulations/predictions)

# Notes
- Works like Pearson's correlation coefficient (r) except Lin's CCC has a
  slope of 1 and intercept of 0. Unlike r which measures only covariation,
  Lin's CCC measures both covariation and correspondence (agreement).
- Lin, L. I. -K. (2000). A note on the concordance correlation coefficient.
  Biometrics, 56, 324–325
- Lin, L. I. -K. (1989). A concordance correlation coefficient to evaluate
  reproducibility. Biometrics, 45(1), 255–268.

# Returns
- `Float64`: Lin's concordance correlation coefficient
"""
function LINCCC(obs::AbstractVector, est::AbstractVector)
    2 * cov(obs, est) / (var(obs) + var(est) + (mean(obs) - mean(est))^2)
end


"""
    AIC(k::Int, order2::Bool=true)

Akaike’s Information Criterion (AIC).

By itself, AIC has no meaning. AIC is meant to be used to compare between models,
where the best model is one with the lowest AIC value.

# Arguments
- `k::Int`: no. of model parameters plus one (see Notes)
- `order2::Bool=true`: `true` for second-order, else `false` for first-order AIC

# Notes
- A simple linear regression equation, y = mx + c, has 3 parameters
  (m and c parameters + 1), so set argument `k` to 3.
- Return a function object that accepts two arguments:
    - `obs::AbstractVector`: vector of observations
    -` est::AbstractVector`: vector of estimations (simulations/predictions)
- Example of how to use:
    a = AIC(3, true)   # 3-parameter equation and 2nd order AIC
    a(obs, est)  # pass the observed and estimated values of the 3-parameter eqn.
- Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference: understanding
  AIC and BIC in Model Selection. Sociological Methods & Research, 33, 261–304.
- Burnham, K. P., & Anderson, D. R. (2002). Model Selection and Multimodel
  Inference: A practical information-theoretic approach (2nd ed.). Springer-Verlag, NY.

# Returns
- `Float64`: Akaike’s Information Criterion (AIC)
"""
function AIC(k::Int, order2::Bool=true)
    function aic_fn(obs::AbstractVector, est::AbstractVector)
        rss = sum((est .- obs) .^ 2)
        n = length(obs)
        mle = rss / n
        -2log(mle) + 2k + (order2 ? 2k * (k + 1) / (n - k - 1) : 0)
    end
end


"""
    BIC(k::Int)

Bayesian information criterion (BIC).

By itself, BIC has no meaning. BIC is meant to be used to compare between models,
where the best model is one with the lowest BIC value.

# Arguments
- `k::Int`: no. of model parameters plus one (see Notes)

# Notes
- A simple linear regression equation, y = mx + c, has 3 parameters
  (m and c parameters + 1), so set argument `k` to 3.
- Return a function object that accepts two arguments:
    - `obs::AbstractVector`: vector of observations
    -` est::AbstractVector`: vector of estimations (simulations/predictions)
- Example of how to use:
    b = BIC(3)   # 3-parameter equation
    b(obs, est)  # pass the observed and estimated values of the 3-parameter eqn.
- Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference: understanding
  AIC and BIC in Model Selection. Sociological Methods & Research, 33, 261–304.
- Burnham, K. P., & Anderson, D. R. (2002). Model Selection and Multimodel
  Inference: A practical information-theoretic approach (2nd ed.). Springer-Verlag, NY.

# Returns
- `Float64`: Bayesian information criterion (BIC)
"""
function BIC(k::Int)
    function bic_fn(obs::AbstractVector, est::AbstractVector)
        rss = sum((est .- obs) .^ 2)
        n = length(obs)
        mle = rss / n
        -2log(mle) + k * log(n)
    end
end

end     # module
