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
    using PlutoUI
    using HDF5, FileIO, Printf, LaTeXStrings
    using Nemo
    using BigCombinatorics
    using Revise
    using TI_quantum

    # initialize the table for factorials
    Factorial_1()

    NMAX = 2
    PATH_FIGS, PATH_DATA = path()
end


begin
    # g2 function and field for different ϕ and r
    r_list_1 = range(0.1, 0.999, NMAX)
    phi_list = range(0, 2*pi, 2)
    alpha_list = range(0.5, 2.0, 2)

    g2_k = zeros(Float64, length(alpha_list), length(phi_list), length(r_list_1))
    g2_mk = zeros(Float64, length(alpha_list), length(phi_list), length(r_list_1))
    g2_cross = zeros(Float64, length(alpha_list), length(phi_list), length(r_list_1))

    n_k = zeros(Float64, length(alpha_list), length(phi_list), length(r_list_1))
    n_mk = zeros(Float64, length(alpha_list), length(phi_list), length(r_list_1))
    n_k_cross = zeros(Float64, length(alpha_list), length(phi_list), length(r_list_1))
    n_mk_cross = zeros(Float64, length(alpha_list), length(phi_list), length(r_list_1))

    a_k_av = zeros(ComplexF64, length(alpha_list), length(phi_list), length(r_list_1))
    a_mk_av = zeros(ComplexF64, length(alpha_list), length(phi_list), length(r_list_1))

    for k in eachindex(alpha_list), i in eachindex(phi_list), j in eachindex(r_list_1)
        ϕ = phi_list[i]
        r = r_list_1[j]

        if r >= 0.4
            num = 20
        elseif r >= 0.3 && r < 0.4
            num = 30
        else
            num = 60
        end

        alpha = alpha_list[k]
        beta = abs((r - 1) / (r + 1)) * exp(im*ϕ) * alpha

        g2_k[k, i, j], n_k[k, i, j] = real.(g2_function(alpha, beta, r, num, "k"))
        g2_mk[k, i, j], n_mk[k, i, j] = real.(g2_function(alpha, beta, r, num, "-k"))
        g2_cross[k, i, j], n_mk_cross[k, i, j], n_k_cross[k, i, j] = real.(g2_function(alpha, beta, r, num, "k,-k"))

        a_k_av[k, i, j] = op_average(alpha, beta, r, num, "a_k")
        a_mk_av[k, i, j] = op_average(alpha, beta, r, num, "a_-k")
    end
end

# Writing data into file

let
    data_dict = Dict("r" => collect(r_list_1), "phi" => collect(phi_list), 
                    "alpha" => collect(alpha_list), 
                    "beta" => alpha_list*abs.((r_list_1.-1)./(r_list_1.+1))',
                    "g2_k" => g2_k, "g2_mk" => g2_mk,
                    "g2_cross" => g2_cross,
                    "n_k" => n_k, "n_mk" => n_mk,
                    "a_k_av_re" => real(a_k_av), "a_k_av_im" => imag(a_k_av),
                    "a_mk_av_re" => real(a_mk_av), "a_mk_av_im" => imag(a_mk_av),
    )

    NAME_PART = "alpha"*"_beta_r_dependent_more_points"*".h5"
    save(PATH_DATA*"g2_field_erasing_"*NAME_PART, data_dict)

    data_dict_loaded = load(PATH_DATA*"g2_field_erasing_"*NAME_PART)
    data_dict_loaded["g2_cross"] == data_dict["g2_cross"]
end
