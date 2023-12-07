begin
    if pwd()[end-9:end] == "TI_quantum"
        PATH_ENV = "."
    else
        PATH_ENV = "../"
    end
    using Pkg
    Pkg.activate(PATH_ENV)
end

begin
	using Plots
	using BenchmarkTools
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

function coherent_state(alpha::Complex, beta::Complex, r::Real,
                        num::Int; precision=4096)
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

function fock_state(n_1::Integer, n_2::Integer, r::Real,
                    num::Int; precision=4096)
    b = FockBasis(2*num)
    (
    sum([
    convert(ComplexF64, c_l_hyper(n_1, n_2, s - min(n_1, n_2), r; precision)) * 
    tensor(fockstate(b, n_1 + s - min(n_1, n_2)),fockstate(b, n_2 + s - min(n_1, n_2)))
    for s = 0:num])
    )
end

begin
    n_terms = 20
	alpha_list = [0.5]# range(0.5, 3.0, 6)
	r_list = range(0.1, 0.999, 20)
	phi_list = range(0.0, pi, 9)
	entropy_ent_vn = zeros(length(alpha_list), length(phi_list), length(r_list))
end

begin
	for i in eachindex(alpha_list), k in eachindex(phi_list), j in eachindex(r_list)
		alpha = alpha_list[i] + 0im
		ϕ = phi_list[k]
		r = r_list[j]
		
		# compute coherent state
		ψ = coherent_state(alpha*exp(1.0im*ϕ), alpha, r, n_terms)
		
		# find the reduced density matrix
		ρ_red = ptrace(tensor(ψ, dagger(ψ)), 1)
		
		# compute entanglement entropy
		# entropy_ent_vn[i, k, j] = real(entropy_vn(ρ_red))
		entropy_ent_vn[i, k, j] = real(entropy_vn(tensor(ψ, dagger(ψ))))
	end
end


entropy_ent_vn[1, :,:]'

let
    gr()
	ind = 1
	f = plot(
		plot(r_list, entropy_ent_vn[1,:,:]'.+10*abs(minimum(entropy_ent_vn)), 
				label=reshape(L"\varphi/\pi = ".*string.(phi_list/pi), 1, length(phi_list)),
				xlabel=L"r = n_0/n_1",
				ylabel=L"S_\mathrm{ent}",
				#title=L"\alpha = 0.5",
                yscale=:log10,
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

let 
    psi_test = coherent_state(1.0*exp(1.0im*pi), 1.0+0im, 0.1, 20)
    rho_test_1 = ptrace(tensor(psi_test, dagger(psi_test)), 1)
    rho_test_2 = ptrace(tensor(psi_test, dagger(psi_test)), 2)
    xvec = range(-5.0, 5.0, 100)
    w_test_1 = wigner(rho_test_1, xvec, xvec)
    w_test_2 = wigner(rho_test_2, xvec, xvec)

    l = @layout([a b])
    p_1 = heatmap(xvec, xvec, w_test_1,
        aspect_ratio = 1
    )
    p_2 = heatmap(xvec, xvec, w_test_2,
        aspect_ratio = 1
    )
    plot(p_1, p_2, layout=l, size=(750,300))
end



# Save data
let
    data_dict = Dict("r" => collect(r_list), "n_terms" => n_terms,
                    "phi" => collect(phi_list),
                    "alpha" => collect(alpha_list), 
                    "beta" => collect(alpha_list),
                    "entropy_ent_vn" => entropy_ent_vn,
    )

    NAME_PART = "nterms"*string(n_terms)*"_equal_alpha_beta"*".h5"
    save(PATH_DATA*"entropy_ent_"*NAME_PART, data_dict)

    data_dict_loaded = load(PATH_DATA*"entropy_ent_"*NAME_PART)
    data_dict_loaded["entropy_ent_vn"] == data_dict["entropy_ent_vn"]
end