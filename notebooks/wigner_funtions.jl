### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 42b2f111-6ae6-4931-908a-3599df2271ce
begin
	using Plots
	using BenchmarkTools
	using PlutoUI
	using HDF5, FileIO, Printf, LaTeXStrings
	using Nemo
	using BigCombinatorics
	using QuantumOptics
	using Revise
	using TI_quantum

	# initialize the table for factorials
	Factorial_1()

	NMAX = 20
	PATH_FIGS, PATH_DATA = path()
end

# ╔═╡ 95cf1d7c-fa49-11ed-0c98-43f3ad8a619b
begin
    if pwd()[end-9:end] == "TI_quantum"
        PATH_ENV = "."
    else
        PATH_ENV = "../"
    end
    using Pkg
    Pkg.activate(PATH_ENV)
end

# ╔═╡ f705723c-d238-4159-b0a2-c00def896a50
begin
	function coherent_state(alpha::Complex, beta::Complex, r::Real, num::Int; precision=4096)
		b = FockBasis(2*num)
		factor_1 = exp(-0.5*(abs(alpha)^2.0 + abs(beta)^2.0))
		factor_2 = (
		sum([convert(ComplexF64, alpha^(n_1) * beta^(n_2) *
	    (1/sqrt(Factorial_1(n_1) * Factorial_1(n_2))) * 
	    convert(ComplexF64, c_l_hyper(n_1, n_2, s - min(n_1, n_2), r; precision))) * 
		tensor(fockstate(b, n_1 + s - min(n_1, n_2)),fockstate(b, n_2 + s - min(n_1, n_2)))
	    for n_1 = 0:num, n_2 = 0:num, s = 0:num])
		)
		factor_1 * factor_2
	end
end

# ╔═╡ 8c9d6292-02cf-41bc-85c4-2a260121f8dd
begin
	function fock_state(n_1::Integer, n_2::Integer, r::Real, num::Int; precision=4096)
		b = FockBasis(2*num)
		(
		sum([
	    convert(ComplexF64, c_l_hyper(n_1, n_2, s - min(n_1, n_2), r; precision)) * 
		tensor(fockstate(b, n_1 + s - min(n_1, n_2)),fockstate(b, n_2 + s - min(n_1, n_2)))
	    for s = 0:num])
		)
	end
end

# ╔═╡ 70808caa-89b1-4aad-b29a-0d1136ff0804
begin
	alpha_list = [0.5]# range(0.5, 3.0, 6)
	r_list = range(0.1, 0.999, 20)
	phi_list = range(0.0, pi, 9)
	entropy_ent_vn = zeros(length(alpha_list), length(phi_list), length(r_list))
end

# ╔═╡ aee2d579-807e-4c04-b494-5117bc7c5956
begin
	for i in eachindex(alpha_list), k in eachindex(phi_list), j in eachindex(r_list)
		alpha = alpha_list[i] + 0im
		ϕ = phi_list[k]
		r = r_list[j]
		
		# compute coherent state
		ψ = coherent_state(alpha*exp(1.0im*ϕ), alpha, r, 30)
		
		# find the reduced density matrix
		ρ_red = ptrace(tensor(ψ, dagger(ψ)), 1)
		
		# compute entanglement entropy
		entropy_ent_vn[i, k, j] = real(entropy_vn(ρ_red))
	end
end

# ╔═╡ b46e9951-add4-4661-84d2-624ec6889ed7
entropy_ent_vn[1, :,:]'

# ╔═╡ 01006c61-c4bf-4531-8d53-83b6560e19f7
let
	ind = 1
	f = plot(
		plot(r_list, entropy_ent_vn[1,:,:]', 
				label=reshape(L"\varphi/\pi = ".*string.(phi_list/pi), 1, length(phi_list)),
				xlabel=L"r = n_0/n_1",
				ylabel=L"S_\mathrm{ent}",
				#title=L"\alpha = 0.5",
				line = 2,
		),
		plot(phi_list/pi, entropy_ent_vn[1,:,ind], 
				label=L"r = ".*string.(round(r_list[ind], digits=2)),
				xlabel=L"\varphi/\pi",
				ylabel=L"S_\mathrm{ent}",
				#title=L"\alpha = 0.5",
				line = 2,
		),
	)
	NAME_PART = "entropy_ent_"
	#savefig(f, PATH_FIGS*NAME_PART*"convergence_equal_alpha_beta.pdf")
	f
end

# ╔═╡ bed891a9-39a2-4004-96ec-fe83443a6440
let
	f = plot(r_list, entropy_ent_vn[1, :], 
				label="α = ".*string.(alpha_list[1]),
				xlabel=L"r = n_0/n_1",
				ylabel=L"S_\mathrm{ent}",
				#title=L"\alpha = 0.5",
				line = 2,
	)
	NAME_PART = "entropy_ent_"
	#savefig(f, PATH_FIGS*NAME_PART*"convergence_equal_alpha_beta_alpha0.5.pdf")
	f
end

# ╔═╡ c692340e-045b-4e5c-8cc1-b4c78527fffa
md"""
## Fock states
"""

# ╔═╡ 75f23d28-8960-4c79-9542-e270ceca4ad1
begin
	n_list = 0:11
	entropy_ent_vn_f = zeros(length(n_list), length(r_list))
end

# ╔═╡ 5c8b0911-2d23-4af4-9875-e46b76873839
begin
	for i in eachindex(n_list), j in eachindex(r_list)
		n_i = n_list[i]
		r = r_list[j]
		
		# compute coherent state
		ψf = fock_state(n_i, n_i, r, 50)
		
		# find the reduced density matrix
		ρf_red = ptrace(tensor(ψf, dagger(ψf)), 1)
		
		# compute entanglement entropy
		entropy_ent_vn_f[i, j] = real(entropy_vn(ρf_red))
	end
end

# ╔═╡ 12ff42e3-2151-4de1-b6df-b5ecc911491f
let
	f = plot(r_list, entropy_ent_vn_f', 
				label=reshape("n = ".*string.(n_list), 1, length(n_list)),
				xlabel=L"r = n_0/n_1",
				ylabel=L"S_\mathrm{ent}",
				#title=L"\alpha = 0.5",
				line = 2,
	)
	NAME_PART = "entropy_ent_"
	#savefig(f, PATH_FIGS*NAME_PART*"convergence_equal_alpha_beta.pdf")
	f
end

# ╔═╡ 0fe7df6c-4903-422f-8362-b55041b67ef6
function entanglement_entropy_f(n_1::Integer, n_2::Integer, r::Real, num::Int; precision=4096)
	entropy = 0.0
	for s = 0:num
		c = c_l_hyper(n_1, n_2, s - min(n_1, n_2), r; precision)
		entropy -= convert(Float64, abs(c)^2*log(abs(c)^2))
	end
	entropy
end

# ╔═╡ 553e21e3-fcef-4764-922b-27a921b64446
let
	r_list_0 = range(1.001, 10.0, 100)
	plot(r_list_0, [entanglement_entropy_f(j, j, r_list_0[i], 200) for i in eachindex(r_list_0), j = 0:11], lw=2,
	xscale=:log10)
	#plot!(r_list, entropy_ent_vn_f[6,:])
end

# ╔═╡ bd39728b-7fcb-4587-b69e-b8bf4e044995
md"""
## Save data
"""

# ╔═╡ 00b7ec0d-95dd-4c1d-8133-6cb69f265080
let
    data_dict = Dict("r" => collect(r_list), "n_terms" => 50,
                    "alpha" => collect(alpha_list), 
                    "beta" => collect(alpha_list),
                    "entropy_ent_vn" => entropy_ent_vn,
    )

    NAME_PART = "nterms50_equal_alpha_beta"*".h5"
    save(PATH_DATA*"entropy_ent_"*NAME_PART, data_dict)

    data_dict_loaded = load(PATH_DATA*"entropy_ent_"*NAME_PART)
    data_dict_loaded["entropy_ent_vn"] == data_dict["entropy_ent_vn"]
end

# ╔═╡ 44282925-839f-4a22-9bf2-472ce820635b
md"""
## Test zone
"""

# ╔═╡ 4ef1f912-fba9-4ccd-8770-b8649da5fe3e
test = coherent_state(0.5+0im, 0.5+0im, 0.147316, 20)

# ╔═╡ 9a829c42-b3f1-4c92-9026-87e4ca768963
test_rho = tensor(test, dagger(test))

# ╔═╡ 7f833634-b866-460b-858b-b16257232e50
test_rho_red = ptrace(test_rho, 1)

# ╔═╡ 6d0c1f65-8c83-4a42-8554-062770a00cdd
# ╠═╡ disabled = true
#=╠═╡
entropy_vn(test_rho)
  ╠═╡ =#

# ╔═╡ d4794a7c-4d97-4bf4-abcf-663ccbd2d301
entropy_vn(test_rho_red)

# ╔═╡ 80870ab3-7ace-4258-9cee-8496b76d9292
aa = coherentstate(FockBasis(16), 2)

# ╔═╡ 0b4a66ec-87ac-4e19-ae27-80d32b03c49e
entropy_vn(tensor(aa, dagger(aa)))

# ╔═╡ 440e6b93-8023-4bd3-b655-eb5fa2d5e3c5
bb = coherentstate(FockBasis(16), 2)

# ╔═╡ 04a97c0f-a2b2-4c5b-b777-7798a81c9712
aabb = tensor(aa, bb)

# ╔═╡ 5540bf96-e67f-411c-af77-0f231b131408
entropy_vn(tensor(aabb, dagger(aabb)))

# ╔═╡ Cell order:
# ╠═95cf1d7c-fa49-11ed-0c98-43f3ad8a619b
# ╠═42b2f111-6ae6-4931-908a-3599df2271ce
# ╠═f705723c-d238-4159-b0a2-c00def896a50
# ╠═8c9d6292-02cf-41bc-85c4-2a260121f8dd
# ╠═70808caa-89b1-4aad-b29a-0d1136ff0804
# ╠═aee2d579-807e-4c04-b494-5117bc7c5956
# ╠═b46e9951-add4-4661-84d2-624ec6889ed7
# ╠═01006c61-c4bf-4531-8d53-83b6560e19f7
# ╠═bed891a9-39a2-4004-96ec-fe83443a6440
# ╟─c692340e-045b-4e5c-8cc1-b4c78527fffa
# ╠═75f23d28-8960-4c79-9542-e270ceca4ad1
# ╠═5c8b0911-2d23-4af4-9875-e46b76873839
# ╠═12ff42e3-2151-4de1-b6df-b5ecc911491f
# ╠═0fe7df6c-4903-422f-8362-b55041b67ef6
# ╠═553e21e3-fcef-4764-922b-27a921b64446
# ╟─bd39728b-7fcb-4587-b69e-b8bf4e044995
# ╠═00b7ec0d-95dd-4c1d-8133-6cb69f265080
# ╟─44282925-839f-4a22-9bf2-472ce820635b
# ╠═4ef1f912-fba9-4ccd-8770-b8649da5fe3e
# ╠═9a829c42-b3f1-4c92-9026-87e4ca768963
# ╠═7f833634-b866-460b-858b-b16257232e50
# ╠═6d0c1f65-8c83-4a42-8554-062770a00cdd
# ╠═d4794a7c-4d97-4bf4-abcf-663ccbd2d301
# ╠═80870ab3-7ace-4258-9cee-8496b76d9292
# ╠═0b4a66ec-87ac-4e19-ae27-80d32b03c49e
# ╠═440e6b93-8023-4bd3-b655-eb5fa2d5e3c5
# ╠═04a97c0f-a2b2-4c5b-b777-7798a81c9712
# ╠═5540bf96-e67f-411c-af77-0f231b131408
