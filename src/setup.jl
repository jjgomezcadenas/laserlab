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
	μW, mW, W


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
- `r::Unitful.Time`   : laser period 
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
end


"""
GaussianLaser

Representation of a Gaussian laser defined by a laser and a waist  

# Fields
- `laser::Laser`       : A Laser 
- `w0::Unitful.Length` : Waist of laser 

"""
struct GaussianLaser 
	laser::Laser
	w0::Unitful.Length
end


"""
	struct Objective

Simple representation of a microscope objective

# Fields
- `name::String`       : identifies the objective
- `NA::Float64`        : Numerical aperture
- `M::Float64`         : Magnification
- `f::typeof(1.0mm)`   : focal distance of the lens
- `d::typeof(1.0mm)`   : diameter of the lens

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


# GaussianLaser functions
"""
	w0_large_zr(laser::Laser, obj::Objective)

returns w0 in the limit of large zr. 
This applies to a (gaussian) laser filling the back lens of a microscope objective.  
 
# Fields

- `laser::Laser` : A laser, filling the back lens of a microscope objective
- `obj::Objective` : The microscope objective  

"""
w0_large_zr(laser::Laser, obj::Objective) = laser.λ/(π * obj.NA)


"""
	gzr(gl::GaussianLaser)

returns zr =  π * w0^2/λ, this is the zr that corresponds to the beam waist w0.
This applies only to a free Gaussian Laser, not to a Gaussian Laser focused on an 
objective (in this case zr is assumed to be very large) 

"""
gzr(gl::GaussianLaser)  = π * gl.w0^2/gl.laser.λ

"""
	gI0(gl::GaussianLaser)

returns I0 =  2.0*P/(π*w0^2)

"""
gI0(gl::GaussianLaser)  = 2.0 * gl.laser.P/(π * gl.w0^2)

"""
    gρ0(gl::GaussianLaser)

returns number of photons corresponding to I0

"""
gρ0(gl::GaussianLaser)  = n_photons(gl.laser.λ, 2.0 * gl.laser.P)/(π * gl.w0^2)


"""
	gwz(gl::GaussianLaser)

returns w(z) = w0 * Sqrt(1 +(z/zr)^2)
# Fields

"""
function gwz(gl::GaussianLaser)
	function wxz(z::Unitful.Length)
		return gl.w0 * sqrt(1 + (z/gzr(gl)^2))
	end
	return wxz
end


"""
	gI(gl::GaussianLaser)

returns I(r,z) = ρ0 [w0/w(z)]^2 * exp(-2(r/w(z)^2))) (Hz/mm^2)

"""
function gI(gl::GaussianLaser)
	function Ix(r::Unitful.Length, z::Unitful.Length)
		wzz = gwz(gl)
		return gρ0(gl) * (gl.w0/wzz(z))^2 * exp(-2*(r/wzz(z))^2)
	end
	return Ix
end


"""
	focus_lens_at_beam_waist(gl::GaussianLaser)

Returns a new GaussianLaser beam after focusing with a lens located at beam waist.
The new beam waist becomes: w' = (f * λ) / (π * w0), where f is the focal.

"""
function focus_lens_at_beam_waist(gl::GaussianLaser, f::Unitful.Length)
	w0 = f * gl.laser.λ/(π * gl.w0)
	return GaussianLaser(gl.laser, w0)
end

