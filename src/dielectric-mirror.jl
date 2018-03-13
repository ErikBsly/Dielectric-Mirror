# Refractive indices
# SiO2: https://refractiveindex.info/?shelf=main&book=Ta2O5&page=Bright-amorphous
# Ta2O5: https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson

function ReflectivityOfDielectricMirror(λ::Float64, λ₀::Float64, θ_::Float64 = 0.0, N::Int = 20)
	θ = deg2rad(θ_)
	n₀, n₁, n₂, ns = 1.0, 2.1411, 1.4599, 1.0 # Sorrounding(e.g. Ethylene Glyc.), Ta2O5, SiO2, Substrate
	
	# arrays
	n = repeat([n₁, n₂], outer = N)
	d = λ₀ ./ 4n
	k = 2π/λ .* n
	δ = k.*d*cos(θ)
	M = [ [	cos(δₓ) 		1im*sin(δₓ)/nₓ;
			1im*sin(δₓ)*nₓ	cos(δₓ)			] for (δₓ, nₓ) in zip(δ, n)]

	# https://de.wikipedia.org/wiki/Bragg-Spiegel (solve final equation component wise), see also Optics Sheet 10
	v₁, v₂ = [n₀ -1; n₀ 1] * prod(M) * [1, ns]
	R = abs2(v₁/v₂)

	return R 
end

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

λ_ = 250.0 : 0.5 : 700.0
T = [1-ReflectivityOfDielectricMirror(λ*1e-9, 550e-9) for λ in λ_] # T = 1 - R

writedlm("out.dat", zip(λ_, T), "\t")
