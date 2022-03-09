using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
using Interpolations

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
	A, N, mol, mmol, V, L, M

	import PhysicalConstants.CODATA2018: N_A


"""
	struct Fov

Represent a field of view

# Fields
- `d::Unitful.Length`  : diameter of Fov
- `z::Unitful.Length`  : thickness
- `a::Unitful.Area`    : area (computed)
- `v::Unitful.Volume`  : volume (computed)

"""
struct Fov
    d::Unitful.Length
    z::Unitful.Length
	a::Unitful.Area
    v::Unitful.Volume

	function Fov(d,z)
		a = π * (d/2.)^2
		v = a * z
		new(d,z,a,v)
	end
end

"""
	struct Fluorophore

Represent a fluorescent molecule

# Fields
- `expeak::Units (nm)`           : peak excitation
- `empeak::Units (nm)`           : peak emission
- `ϵ::Units(``cm^{-1} M^{-1}``)` : molar extinction coefficient
- `Q::Float64`                   : Quantum efficiency
- `σ::Units(``cm^2``)`           : attenuation cross section
"""
struct Fluorophore
    expeak::typeof(1.0nm)
    empeak::typeof(1.0nm)
    ϵ::typeof(1.0/(cm*M))
    Q::Float64
	σ::typeof(1.0cm^2)
	function Fluorophore(ex, en, ϵ, Q)
		σ = log(10) * uconvert(cm^2/mol, ϵ) / N_A
		new(ex, en, ϵ, Q, σ)
	end
end


"""
	struct Monolayer

Represents a monolayer

# Fields
- `p::Unitful.Length`  : pitch between molecules
- `σ::Units (``molecules/μm^2``)`  : number of molecules per area
"""
struct Monolayer
	p::Unitful.Length 
    σ::typeof(1.0*μm^-2)

	function Monolayer(p)
		σ = uconvert(μm^-2, 1.0/p^2)
		new(p, σ)

	end
end


abstract type Laser end

"""
	struct CLaser <: Laser

Simple representation of a continous laser

# Fields
- `λ::typeof(1.0nm)`  : Laser wavelength
- `P::typeof(1.0mW)`  : Power

"""
struct CLaser <: Laser
	λ::Unitful.Length
	P::Unitful.Power
end


"""
	struct PulsedLaser

Simple representation of a pulsed laser

# Fields
- `λ::Unitful.Length`  : Laser wavelength
- `Pk::Unitful.Power`  : peak Power
- `P::Unitful.Power`  : average power
- `f::Unitful.Frequency`  : laser frequency
- `w::Unitful.Time` : pulse width 

"""
struct PulsedLaser  <: Laser
	λ::Unitful.Length
	Pk::Unitful.Power
	f::Unitful.Frequency
	P::Unitful.Power
	w::Unitful.Time
	

	function PulsedLaser(λ::Unitful.Length,
						 Pk::Unitful.Power,
						 f::Unitful.Frequency,
		                 w::Unitful.Time)

		
		P = uconvert(mW, Pk * w * f)
		new(λ, Pk, f, P,  w)
	end

	function PulsedLaser(λ::Unitful.Length,
						 Pk::Unitful.Power,
						 f::Unitful.Frequency,
						 P::Unitful.Power,
						 w::Unitful.Time)


			new(λ, Pk, f, P,  w)
	end
end


# """
# GaussianLaser

# Representation of a Gaussian laser defined by a laser and a waist  

# # Fields
# - `laser::Laser`       : A Laser 
# - `w0::Unitful.Length` : Waist of laser 

# """
# struct GaussianLaser 
# 	laser::Laser
# 	w0::Unitful.Length
# end


"""
	struct TLens

Represents a thin lens

# Fields
- `M::Real`         : Magnification
- `f::Unitful.Length`   : focal distance of the lens
- `d::Unitful.Length`   : diameter of the lens

"""
struct TLens
    M::Real 
	f::Unitful.Length
	d::Unitful.Length
	
end

"""
	propagate(tl::TLens, y1::Real, z1::Real)

Propagate ray from (z1, y1) to (z2, y2) through thin lens tl 

# Fields
- `tl::TLens`         : A thin lens 
- `z1::Unitful.Length`   : point z1 
- `y1::Unitful.Length`   : point y1 

"""
function propagate(tl::TLens, y1::Unitful.Length, z1::Unitful.Length)
	y2 = tl.M * y1
	iz2 = (1.0/tl.f) - (1.0/z1)
	y2, 1.0/iz2
end


"""
	struct Objective

Simple representation of a microscope objective

# Fields
- `name::String`       : identifies the objective
- `NA::Float64`        : Numerical aperture
- `M::Float64`         : Magnification
- `f::typeof(1.0mm)`   : focal length 
- `d::typeof(1.0mm)`   : entrance pupil diameter

"""
struct Objective
	name::String
    NA::Real 
    M::Real 
	f::Unitful.Length
	d::Unitful.Length

	function Objective(name::String, NA::Float64,  M::Float64)
		new(name, NA, M, -1.0mm, -1.0mm)
	end

	function Objective(name::String, f::Unitful.Length, d::Unitful.Length, M::Float64)
		ff = f / d   # https://www.eckop.com/resources/optics/numerical-aperture-and-f-number/
		NA = 1.0/(2.0*ff)
		new(name, NA, M, f,d)
	end
end



#CCD is defined in terms of a function which returns the CCD response
#(e.g, function ccd returns a response function, which gives efficiency
# as a function of wavelength)
"""
	ccd(lmin::Float64=350.0, lmax::Float64=1000.0)

Return the efficiency of a generic CCD as a function of wavelength.

# Fields

- `lmin::Float64=350.0` : Minimum wavelength for which efficiency is defined
- `lmax::Float64=350.0` : Maximum wavelength for which efficiency is defined

"""
function ccd(lmin::Float64=350.0, lmax::Float64=1000.0)
	function eff(l::Float64)
		if l < lmin || l > lmax
			return 0.
		else
			wl = 350.:50.:1000.
			ϵ = [0.3, 0.4,0.65,0.78,0.82,0.82,0.8,0.72,0.62,0.5,0.37,
			  0.24,0.12,0.07]
			e = CubicSplineInterpolation(wl, ϵ)
			return e(l)
		end
	end
	return eff
end


#FUNCTIONS

"""
	photon_energy(λ::Unitful.Length)

Given wavelength of photon return its energy.
# Fields

- `λ::Unitful.Length`  : Photon wavelength

"""
function photon_energy(λ::Unitful.Length)
	uconvert(eV, λ, Spectral())
end


"""
	delivered_energy(laser::Laser, t::Unitful.Time)

Delivered energy of a laser in a given time.
# Fields

- `laser::Laser`     : Laser
- `t::Unitful.Time`  : Time in which target is illuminated

"""
function delivered_energy(laser::Laser, t::Unitful.Time)
	laser.P * t
end

"""
	n_photons(laser::Laser)

Rate of photons (number of photons per unit time) produced by a laser
# Fields

- `laser::Laser`     : Laser
- `t::Unitful.Time`  : Time in which target is illuminated

"""
function n_photons(laser::Laser)
	uconvert(Hz, laser.P / photon_energy(laser.λ))
end


"""
	n_photons(λ::Unitful.Length, p::Unitful.Power)

Rate of photons (number of photons per unit time) corresponding to a wavelength
λ and a power P

# Fields

- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power

"""
function n_photons(λ::Unitful.Length, p::Unitful.Power)
	uconvert(Hz, p / photon_energy(λ))
end


"""
	n_photons_int(laser::Laser, t::Unitful.Time)

Integrated number of photons in a given time emitted by a laser

# Fields

- `laser::Laser`    : Laser
- `t::Unitful.Time` : time of measurement

"""
function n_photons_int(laser::Laser, t::Unitful.Time)
	uconvert(eV,delivered_energy(laser, t)) / photon_energy(laser.λ)
end


"""
	photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area)

number of photons per unit time per unit area

# Fields

- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power
- `a::Unitful.Area`   : Area

"""
function photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area)
	return n_photons(λ, p)/ a
end


"""
	photon_density(l::Laser, fov::Fov)

number of photons per unit time per unit area, in a Fov illuminated by a laser

# Fields

- `laser::Laser` : Laser
- `fov::Fov`     : Field of view

"""
function photon_density(laser::Laser, fov::Fov)
	return n_photons(laser) / fov.a
end


"""
	diffraction_limit(l::Laser, obj:: Objective)

Return the diameter of diffractive spot for a laser l
focused with an objective obj

# Fields

- `l::Laser`         : A laser
- `obj:: Objective`  : An objective
"""
function diffraction_limit(l::Laser, obj:: Objective)
    return 1.83 * l.λ/(2 * obj.NA)
end


"""
	geometrical_acceptance(d::Float64, D::Float64)

Compute the fraction of photons that make it through an iris
of diameter D located at a distance d from the emission point.

# Fields

- `d::Float64`   : distance between emission point and iris
- `D::Float64`   : Diameter of iris

"""
function geometrical_acceptance(d::Float64, D::Float64)
	return 0.5(1. - d/sqrt(d^2 + (D/2.)^2))
end


"""
	transmission(objective::Objective)

Compute the transmission of an objective (depends only of NA).

# Fields

- `objective::Objective` : Objective

"""
function transmission(objective::Objective)
	A = objective.NA
	A <=1 ? (1 - sqrt(1 - A^2)) /2 : 0.5
end


"""
	transmission(A::Float64)

Compute the transmission as a function of NA.

# Fields

- `A::Float64` : Numerical acceptance (NA)

"""
function transmission(A::Float64)
	A <=1 ? (1 - sqrt(1 - A^2)) /2 : 0.5
end


# fluorophores

"""
	fluorescence(f::Fluorophore, I::Quantity)

Number of  photons emitted per unit time when fluorosphore F
is illuminated with a laser of photon density I

# Fields

- `f::Fluorophore`    : Define fluorophore
- `I::(``nγ/cm^2``)`  : Photon density
"""
fluorescence(f::Fluorophore, I::typeof(1.0Hz/cm^2)) =  f.σ * f.Q * I
	