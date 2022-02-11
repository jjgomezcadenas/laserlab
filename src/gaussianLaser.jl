using Unitful


import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W


    """
	GaussianLaserObj

Representation of a Gaussian laser that fills the objective lens 
assuming zR >> f, where f is the focal 

# Fields
- `laser::Laser`       : A laser type
- `obj::Objective`     : An objective type
- `w0::typeof(1.0nm)`  : Waist of laser at focusing point
- `I0::typeof(1.0mW/cm^2)` : Intensity at r = 0, z= 0
- `ρ0::typeof(1.0Hz/cm^2)` : Photon density at r = 0, z= 0


"""
struct GaussianLaserObj
	laser::Laser
	obj::Objective
	w0::Unitful.Length
	I0::typeof(1.0mW/cm^2)
	ρ0::typeof(1.0Hz/cm^2)

	function GaussianLaserObj(laser,obj)
		w0 = laser.λ/(π * obj.NA)
		I0 = uconvert(mW/cm^2, 2.0 * laser.P/(π*w0^2))
		ρ0 = uconvert(Hz/cm^2, n_photons(laser)/(π*w0^2))
		new(laser,obj, w0, I0, ρ0)
	end
end

"""
	GLaser

Representation of a Gaussian laser defined by wavelength, power and wait 

# Fields
- `λ::typeof(1.0nm)`       : Laser wavelength
- `P::typeof(1.0mW)`       : Power
- `w0::typeof(1.0nm)`      : Waist of laser at focusing point
- `zr::typeof(1.0nm)`      : Rayleigh length
- `I0::typeof(1.0mW/cm^2)` : Intensity at r = 0, z= 0
- `ρ0::typeof(1.0Hz/cm^2)` : Photon density at r = 0, z= 0


"""

function glaser_(λ::Unitful.Length, P::Unitful.Power, w0::Unitful.Length)
	zr =  π * w0^2/λ
	I0 =  2.0*P/(π*w0^2)
	ρ0 =  n_photons(λ, 2.0*P)/(π*w0^2)
	return zr, I0, ρ0
end

mutable struct GLaser
	λ::Unitful.Length
	P::Unitful.Power
	w0::Unitful.Length
end

# GLaser functions

"""
	gzr(gl::GLaser)

returns zr =  π * w0^2/λ

"""
gzr(gl::Glaser)  = π * gl.w0^2/gl.λ

"""
	gI0(gl::GLaser)

returns I0 =  2.0*P/(π*w0^2)

"""
gI0(gl::Glaser)  = 2.0 * gl.P/(π * gl.w0^2)

"""
    gρ0(gl::GLaser)

returns number of photons corresponding to I0

"""
gρ0(gl::Glaser)  = n_photons(gl.λ, 2.0 * gl.P)/(πg * gl.w0^2)


"""
	gwz(gl::GLaser)

returns w(z) = w0 * Sqrt(1 +(z/zr)^2)
# Fields

"""
function gwz(gl::GLaser)
	function wxz(z::Unitful.Length)
		return gl.w0 * sqrt(1 + (z/gl.zr)^2)
	end
	return wxz
end


"""
	gI(gl::GLaser)

returns I(r,z) = ρ0 [w0/w(z)]^2 * exp(-2(r/w(z)^2))) (Hz/mm^2)

"""
function gI(gl::GLaser)
	function Ix(r::Unitful.Length, z::Unitful.Length)
		wzz = gwz(gl)
		return gl.ρ0 * (gl.w0/wzz(z))^2 * exp(-2*(r/wzz(z))^2)
	end
	return Ix
end

"""
	focus_lens_at_beam_waist(gl::GLaser)

Returns a new GLaser beam after focusing with a lens located at beam waist.
The new beam waist becomes: w' = (f * λ) / (π * w0), where f is the focal.

"""
function focus_lens_at_beam_waist(gl::GLaser, f::typeof(1.0mm))
	w0 = f*gl.λ/(π*gl.w0)
	return GLaser(gl.λ, gl.P, w0)
end
