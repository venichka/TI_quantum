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
    using BenchmarkTools, ProgressMeter
    using PlutoUI
    using HDF5, FileIO, Printf, LaTeXStrings
    using Nemo
    using BigCombinatorics
    using Revise
    using TI_quantum

    # initialize the table for factorials
    Factorial_1()

    PATH_FIGS, PATH_DATA = path()
end

begin
    NMAX = 50
    NMAX_r = 30
    prob_ab_NMAX = zeros(NMAX, NMAX, NMAX_r)
    r_list = range(0.01, 0.99, NMAX_r)
    # r_list = 10.0.^(range(0.005, 2, NMAX_r))
    # r_list = range(1.01, 100, NMAX_r)
    alpha = 5.0
    beta = 5.0
    phase_a = [0, 1]  # [0, 1] -- exp(i π * 0 / 1)
    num = 100
    print("done")
end

begin
    lk = ReentrantLock()
    p = Progress(NMAX*NMAX*NMAX_r)
    Threads.@threads for iijjkk = CartesianIndices((NMAX, NMAX, NMAX_r))
        ii, jj, kk = Tuple(iijjkk)[1], Tuple(iijjkk)[2], Tuple(iijjkk)[3]
        n_1 = ii - 1
        n_2 = jj - 1
        r = r_list[kk]
        lock(lk) do
            result = p_ab(alpha, beta, n_1, n_2, r, num; 
                          phase_fraction_a=phase_a)
            prob_ab_NMAX[ii, jj, kk] = result
            next!(p)
        end
    end
    finish!(p)
end

sum(prob_ab_NMAX[:,:,3])

# Check: mean number of photons in the 2d mode
sum([sum(prob_ab_NMAX[:,:,19], dims=2)[i]*(i-1) for i = 1:NMAX])
plot(r_list, 
    [sum([sum(prob_ab_NMAX[:,:,k], dims=2)[i]*(i-1) for i = 1:NMAX])
    for k in eachindex(r_list)], 
    # xscale=:log10,
    )
plot(r_list, 
    [sum(prob_ab_NMAX[:,:,k])
    for k in eachindex(r_list)], 
    # xscale=:log10,
    )

let 
    idx = 30
	
	xs = [string(i) for i in 0:NMAX-1]
	ys = xs
    z = prob_ab_NMAX[:,:,idx]

    heatmap(xs, ys, z, aspect_ratio = :equal, 
			xlabel = L"N_{k}",
			ylabel = L"N_{-k}",
			framestyle=:box,
			colorbar=:right,
			colorbar_formatter=:scientific,
			title= L"P(N_k, N_{-k})",
			titlefontsize=20,
			fontfamily="serif",
			guidefontsize=18,
			tickfontsize=12,
			colorbar_tickfontsize=12,
			dpi=300,
			rightmargin=10Plots.mm,
  			leftmargin=5Plots.mm,
		)
		annotate!([
		#(xs[35], ys[107],
		#(xs[6], ys[26],
		(xs[4], ys[17],
		Plots.text(L"r = "*string(round(r_list[idx], digits=2)), 
			14, :white, :center))])
end

# Writing into a file

let
    data_dict = Dict("r" => collect(r_list), "n_1" => collect(0:NMAX-1),
                    "n_2" => collect(0:NMAX-1),
                    "alpha" => real(alpha),
                    "beta" => real(beta),
                    "num_sum" => num,
                    "p_ab" => prob_ab_NMAX,
    )

    NAME_PART = "alpha"*string(real(alpha))*"_beta"*string(real(beta))*"_numsum"*string(num)*"_more_points.h5"
    # save(PATH_DATA*"coh_photon_distrib_r_larger1_"*NAME_PART, data_dict)
    save(PATH_DATA*"coh_photon_distrib_"*NAME_PART, data_dict)

    # data_dict_loaded = load(PATH_DATA*"coh_photon_distrib_r_larger1_"*NAME_PART)
    data_dict_loaded = load(PATH_DATA*"coh_photon_distrib_"*NAME_PART)
    data_dict_loaded["p_ab"] == data_dict["p_ab"]
end


