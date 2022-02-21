using Unitful


import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W


	"""
	GaussianLaser 
	
	Representation of a Gaussian laser defined by a laser, the location of the waist (z0) and the waist radius (w0)  
	
	# Fields
	- `laser::Laser`       : A Laser  
	- `w0::Unitful.Length` : radius of the waist  
	- `z0::Unitful.Length` : location of the waist (computed from w0)
	- `I0::typeof(1.0mW/cm^2)` : Intensity of the beam at (0,0) (computed from w0)
	- `θ0::Real ` : divergence of the beam  (computed from w0)
	
	"""
	struct GaussianLaser 
		laser::Laser
		w0::Unitful.Length
		z0::Unitful.Length
		I0::typeof(1.0mW/cm^2)
		γ0::typeof(1.0Hz/cm^2)
		θ0::Real 
		
		function GaussianLaser(laser, w0)
			z0 = π * w0^2 / laser.λ
			I0 = 2 * laser.P / (π * w0^2)
			nγ = uconvert(Hz, laser.P / photon_energy(laser.λ))
			γ0 = 2 * nγ / (π * w0^2)
			θ0 = w0/z0
			new(laser, w0, z0, I0, γ0, θ0)
		end
	end



	"""
	W(gl::GaussianLaser)

returns the beam width: ``W(z) = W_0 \\sqrt{1 + (z/z_0)^2}``

# Fields
	- `gl::GaussianLaser`       : A gaussian Laser 

"""
function W(gl::GaussianLaser)
	Wz(z::Unitful.Length) = gl.w0 * sqrt(1.0 + (z/gl.z0)^2)
	return Wz
end


"""
R(gl::GaussianLaser)

returns the beam radius of curvature: ``R(z) = z ( 1 + (z_0/z)^2) ``

# Fields
- `gl::GaussianLaser`       : A gaussian Laser 

"""
function R(gl::GaussianLaser)
	Rz(z::Unitful.Length) = z * (1.0 + (gl.z0/z)^2)
	return Rz
end


"""
I(gl::GaussianLaser)

returns the beam Intensity: ``I(\\rho, z) = I_0 ( W_0 / W(z))^2 \\exp{-2 \\rho^2/W^2(z)}  ``

# Fields
- `gl::GaussianLaser`       : A gaussian Laser 

"""
function I(gl::GaussianLaser)
	Wz = W(gl)
	Irz(ρ::Unitful.Length, z::Unitful.Length) = gl.I0 * (gl.w0 / Wz(z))^2 * exp(-2.0 * ρ^2/Wz(z)^2)
	return Irz
end


"""
	spot_size(gl::GaussianLaser)

returns 2*w0 that is the diameter of the beam waist, also called spot size   
 
# Fields

- `gl::GaussianLaser` : A gaussian laser 

"""
spot_size(gl::GaussianLaser)  = 2 * gl.w0


"""
	angular_divergence(gl::GaussianLaser)

returns 2θ that is the beam angular divergence    
 
# Fields

- `gl::GaussianLaser` : A gaussian laser 

"""
angular_divergence(gl::GaussianLaser)  = 2 * gl.θ0


"""
	depth_of_focus(gl::GaussianLaser)

returns 2z0   
 
# Fields

- `gl::GaussianLaser` : A gaussian laser 

"""
depth_of_focus(gl::GaussianLaser)  = 2 * gl.z0


"""
	propagate_paralel_beam(laser::Laser, obj::Objective)

returns a gaussian laser with parameter w0 resulting from focusing the paralell laser beam 
filling the entrance pupil of the objective.    
 
# Fields

- `laser::Laser` : A  laser filling the entrance pupil of the objective (assuming parallell rays)
- `obj::Objective` : An objective, used to focuse the laser in the sample. 

"""
function propagate_paralell_beam(laser::Laser, obj::Objective) 
	w0 = laser.λ / (π * obj.NA)
	GaussianLaser(laser, w0)
end

