using LsqFit
using Statistics
using StatsBase

"""
fit_straight_line(tdata::Vector{Float64}, vdata::Vector{Float64}; pa0=[0.0, 0.5], i0=1)

# Fields 

`tdata`: time-like (x) data
`vdata`` : voltage-like (y) data 
`pa0`: initial values of parameters 
`i0` : index of first point in array to be used in the fit 

# Returns 
- coefficient of the field (a vector)
- standard errors of coefficients (a vector)
- yfit = fit(tdata), that is the prediction of the values corresponding to each tdata
  for the fit parameters 

"""
function fit_straight_line(tdata::Vector{Float64}, vdata::Vector{Float64}; pa0=[0.0, 0.5], i0=1)
	tfun(t, a, b) = a + b * t
	pfun(t, p) = p[1] .+ p[2] .* t
	il = length(tdata)
	fit = curve_fit(pfun, tdata[i0:il], vdata[i0:il], pa0)
	coef(fit), stderror(fit), tfun.(tdata, coef(fit)...)
end


"""
fit_expo(tdata::Vector{Float64}, vdata::Vector{Float64}; pa0=[0.0, 0.5], i0=1)

# Fields 

`tdata`: time-like (x) data
`vdata`` : voltage-like (y) data 
`pa0`: initial values of parameters 
`i0` : index of first point in array to be used in the fit 

# Returns 
- coefficient of the field (a vector)
- standard errors of coefficients (a vector)
- yfit = fit(tdata), that is the prediction of the values corresponding to each tdata
  for the fit parameters 

"""
function fit_expo(tdata::Vector{Float64}, vdata::Vector{Float64}; pa0=[0.0, 0.5], i0=1)
	tfun(t, N, λ) = N*exp(-t/λ)
	pfun(t, p) = p[1] * exp.(-t/p[2])
	il = length(tdata)
	fit = curve_fit(pfun, tdata[i0:il], vdata[i0:il], pa0)
	coef(fit), stderror(fit), tfun.(tdata, coef(fit)...)
end


"""
    fit_time_histo(htime::Histogram; pa0=[100.0, 0.5], i0=1)
    
# Fields

	`htime` : An histogram of times 
	`pa0`   : initial value of fit parameters 
	`i0`    : first bin to fit  

# Returns

# Returns 
- coefficient of the field (a vector)
- standard errors of coefficients (a vector)
- yfit = fit(tdata), that is the prediction of the values corresponding to each tdata
  for the fit parameters 

"""
function fit_time_histo(htime::Histogram; pa0=[100.0, 0.5], i0=1)
	tdata = centers(htime)
	vdata = htime.weights
    fit_expo(tdata, vdata; pa0=pa0, i0=i0)
end