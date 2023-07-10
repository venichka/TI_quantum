### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fef9cfab-5ce3-4bd4-8f2f-fd8ba840e06b
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

	NMAX = 20
	PATH_FIGS, PATH_DATA = path()
end

# ╔═╡ ed82f1be-e061-11ed-0a27-757649bc8226
begin
    if pwd()[end-9:end] == "TI_quantum"
        PATH_ENV = "."
    else
        PATH_ENV = "../"
    end
    using Pkg
    Pkg.activate(PATH_ENV)
end

# ╔═╡ e9703909-861e-4c3e-858f-130f7e08d869
md"""
# $g^{(2)}(0)$ analysis. Coherent states
"""

# ╔═╡ 8e4f0f92-2284-4408-9ea9-4a2085753643
begin
	function variance_n(g2, av_n)
		(g2 .- 1) .* av_n.^2 .+ av_n
	end
end

# ╔═╡ 3b52bf51-8e85-4b2f-bfac-a752ee0658ac
begin
	# Load data
	#NAME_PART_full = "alpha"*"_beta_r_dependent_more_points_n"*string(NMAX)*".h5"
	NAME_PART_full = "alpha"*"_beta_r_dependent_more_points.h5"
	data_dict_g2_full = load(PATH_DATA*"g2_field_equal_phi_rgg1_"*NAME_PART_full)
end

# ╔═╡ 19b4d18a-0014-486d-9809-d2ac1613393d
@bind alpha_slider PlutoUI.Slider(1:3, default=1)

# ╔═╡ 5080ed27-3963-4a70-80c8-a64870189189
@bind phi_slider PlutoUI.Slider(1:9, default=1)

# ╔═╡ bed0ed98-13b5-4e7b-9ffc-f4c52bec3dd0
let
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	phi_list = data_dict_g2_full["phi"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	plot(
		plot(r_list, 
			[
				g2_cross[alpha_slider, phi_slider, :],
				g2_k[alpha_slider, phi_slider, :],
				g2_mk[alpha_slider, phi_slider, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = "*string(alpha_list[alpha_slider]),
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[alpha_slider, phi_slider, :] .+ n_mk[alpha_slider, phi_slider, :],
				n_k[alpha_slider, phi_slider, :],
				n_mk[alpha_slider, phi_slider, :]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			title=L"\varphi = "*string(phi_list[phi_slider]/pi)*L"\pi",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[alpha_slider, phi_slider, :], n_k[alpha_slider, phi_slider, :]),
				variance_n(g2_mk[alpha_slider, phi_slider, :], n_mk[alpha_slider, phi_slider, :])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 0.5",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		leg = :topright,
		#size = (800, 600),
	)
end

# ╔═╡ 92cac3ba-455f-45f3-a883-c753c053fc96
let

	ind = [[1, 9]] # alpha: 1:3 -- 0.5,1.25,2; phi: in the range 1:9 -- 0:π
	
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	phi_list = data_dict_g2_full["phi"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	fig = plot(
		plot(r_list, 
			[
				g2_cross[CartesianIndex.(Tuple.(ind)), :][:],
				g2_k[CartesianIndex.(Tuple.(ind)), :][:],
				g2_mk[CartesianIndex.(Tuple.(ind)), :][:]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			#title=L"\alpha = "*string(alpha_list[idd[1][1]]),
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[CartesianIndex.(Tuple.(ind)), :][:] .+ n_mk[CartesianIndex.(Tuple.(ind)), :][:],
				n_k[CartesianIndex.(Tuple.(ind)), :][:],
				n_mk[CartesianIndex.(Tuple.(ind)), :][:]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = "*string(alpha_list[ind[1][1]])*L"\; ; \;\;\;"*L"\varphi = "*string(phi_list[ind[1][2]]/pi)*L"\pi",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[CartesianIndex.(Tuple.(ind)), :][:], n_k[CartesianIndex.(Tuple.(ind)), :][:]),
				variance_n(g2_mk[CartesianIndex.(Tuple.(ind)), :][:], n_mk[CartesianIndex.(Tuple.(ind)), :][:])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 0.5",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		xscale=:log10,
		#leg = :topright,
		layout = (1, 3),
		size = (700, 220),
		bottom_margin=5.0*Plots.mm,
		left_margin=5.0*Plots.mm,
	)
	
	NAME_PART = "g2_field_equal_erasing_"
	#savefig(fig, PATH_FIGS*NAME_PART*"phi"*string(phi_list[ind[1][2]]/pi)*"pi_alpha"*string(alpha_list[ind[1][1]])*".pdf")
	fig
end

# ╔═╡ 51580c92-f227-4d2a-a378-874ba047294c
md"""
 $r$:  $(@bind r_slider PlutoUI.Slider(1:NMAX, default=NMAX))

 $t$:  $(@bind t_slider PlutoUI.Slider(0:400, default=0))

 $x$:  $(@bind x_slider PlutoUI.Slider(0:400, default=0))

 $\alpha$:  $(@bind vphi_slider PlutoUI.Slider(1:3, default=1))
"""

# ╔═╡ 6cc9021d-2958-4778-acfc-2c189e6454a3
let 
	a_k_av = (data_dict_g2_full["a_k_av_re"] .+ 1.0im*data_dict_g2_full["a_k_av_im"])
	a_mk_av = (data_dict_g2_full["a_mk_av_re"] .+ 1.0im*data_dict_g2_full["a_mk_av_im"])

	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	var_n_k = variance_n(g2_k, n_k)
	var_n_mk = variance_n(g2_mk, n_mk)

	e_k = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[r_slider], a_k_av)
	e_mk = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[r_slider], a_mk_av)
	e_full = TI_quantum.e_field_average(x_slider, 2*pi, t_slider, 2*pi*r_list[r_slider], a_k_av, a_mk_av)

	f = plot(
		plot(phi_list/pi, 
			[
				n_k[1, :, r_slider] .+ n_mk[1, :, r_slider],
				n_k[1, :, r_slider],
				n_mk[1, :, r_slider],
				2*[alpha_list[1]^2 for j in eachindex(phi_list)],
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init"], 
			xlabel=L"\varphi/\pi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 2",
			#yscale=:log10,
			color = [:red :black :green :blue], 
			line = ([2 2 2], [:solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[1, :, r_slider].+var_n_mk[1, :, r_slider] var_n_k[1, :, r_slider] var_n_mk[1, :, r_slider]]),
		),
		leg = :topright,
		#layout = (1, 3),
		#size = (700, 220),
	)
end

# ╔═╡ 68511105-2e94-4bcd-b23a-e2fbc2411f5a
let 
	ind = 3
	
	a_k_av = (data_dict_g2_full["a_k_av_re"] .+ 1.0im*data_dict_g2_full["a_k_av_im"])
	a_mk_av = (data_dict_g2_full["a_mk_av_re"] .+ 1.0im*data_dict_g2_full["a_mk_av_im"])

	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	var_n_k = variance_n(g2_k, n_k)
	var_n_mk = variance_n(g2_mk, n_mk)

	e_k = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av)
	e_mk = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_mk_av)
	e_full = TI_quantum.e_field_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av, a_mk_av)

	omega = r_list[ind]
	

	fig = plot(
		plot(phi_list/pi, 
			[
				n_k[1, :, ind] .+ n_mk[1, :, ind] .+ 0.5,
				n_k[1, :, ind] .+ 0.5,
				n_mk[1, :, ind] .+ 0.5,
				(2*[alpha_list[1]^2 for j in eachindex(phi_list)] .+ 0.5)/omega,
			].*omega,
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init energy"], 
			#yscale=:log2,
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\hbar \omega (\langle n \rangle + 1/2)",
			title=L"\alpha = 2",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[1, :, ind].+var_n_mk[1, :, ind] var_n_k[1, :, ind] var_n_mk[1, :, ind]]),
		),
		#leg = :topright,
		layout = (1, 3),
		size = (400, 220),
		bottom_margin=5.0*Plots.mm,
		left_margin=2.0*Plots.mm,
	)
	NAME_PART = "energy_equal_fields_erasing_rgg1_"
	#savefig(fig, PATH_FIGS*NAME_PART*"r"*string(round(r_list[ind], digits=2))*".pdf")
	fig
end

# ╔═╡ 7989d94e-a548-4619-adc3-321ff80d0858
let 
	ind = 20
	
	a_k_av = (data_dict_g2_full["a_k_av_re"] .+ 1.0im*data_dict_g2_full["a_k_av_im"])
	a_mk_av = (data_dict_g2_full["a_mk_av_re"] .+ 1.0im*data_dict_g2_full["a_mk_av_im"])

	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	var_n_k = variance_n(g2_k, n_k)
	var_n_mk = variance_n(g2_mk, n_mk)

	e_k = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av)
	e_mk = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_mk_av)
	e_full = TI_quantum.e_field_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av, a_mk_av)

	omega = r_list[ind]
	

	fig = plot(
		plot(phi_list/pi, 
			[
				n_k[1, :, ind] .+ n_mk[1, :, ind],
				n_k[1, :, ind],
				n_mk[1, :, ind],
				(2*[alpha_list[1]^2 for j in eachindex(phi_list)]),
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init sum"], 
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 2",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[1, :, ind].+var_n_mk[1, :, ind] var_n_k[1, :, ind] var_n_mk[1, :, ind]]),
		),
		#leg = :topright,
		layout = (1, 3),
		size = (700, 220),
		bottom_margin=5.0*Plots.mm,
		left_margin=2.0*Plots.mm,
	)
	NAME_PART = "n_equal_fields_erasing_rgg1_"
	#savefig(fig, PATH_FIGS*NAME_PART*"r"*string(round(r_list[ind], digits=2))*".pdf")
	fig
end

# ╔═╡ 216f1946-b172-4b0b-b6e2-0f64015f88f9
let 
	a_k_av = (data_dict_g2_full["a_k_av_re"] .+ 1.0im*data_dict_g2_full["a_k_av_im"])
	a_mk_av = (data_dict_g2_full["a_mk_av_re"] .+ 1.0im*data_dict_g2_full["a_mk_av_im"])

	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	var_n_k = variance_n(g2_k, n_k)
	var_n_mk = variance_n(g2_mk, n_mk)

	e_k = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[r_slider], a_k_av)
	e_mk = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[r_slider], a_mk_av)
	e_full = TI_quantum.e_field_average(x_slider, 2*pi, t_slider, 2*pi*r_list[r_slider], a_k_av, a_mk_av)

	f = plot(
		plot(phi_list/pi, 
			[
				n_k[1, :, r_slider] .+ n_mk[1, :, r_slider],
				n_k[1, :, r_slider],
				n_mk[1, :, r_slider],
				2*[alpha_list[1]^2 for j in eachindex(phi_list)],
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init"], 
			xlabel=L"\varphi/\pi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 0.5",
			#yscale=:log10,
			color = [:red :black :green :blue], 
			line = ([2 2 2], [:solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[1, :, r_slider].+var_n_mk[1, :, r_slider] var_n_k[1, :, r_slider] var_n_mk[1, :, r_slider]]),
		),
		plot(phi_list/pi, 
			[
				n_k[2, :, r_slider] .+ n_mk[2, :, r_slider],
				n_k[2, :, r_slider],
				n_mk[2, :, r_slider],
				2*[alpha_list[2]^2 for j in eachindex(phi_list)],
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init"], 
			xlabel=L"\varphi/\pi",
			ylabel=L"\langle n \rangle",
			#yscale=:log10,
			title=L"\alpha = 1.25",
			color = [:red :black :green :blue], 
			line = ([2 2 2], [:solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[2, :, r_slider].+var_n_mk[2, :, r_slider] var_n_k[2, :, r_slider] var_n_mk[2, :, r_slider]]),
		),
		plot(phi_list/pi, 
			[
				n_k[3, :, r_slider] .+ n_mk[3, :, r_slider],
				n_k[3, :, r_slider],
				n_mk[3, :, r_slider],
				2*[alpha_list[2]^2 for j in eachindex(phi_list)],
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init"], 
			xlabel=L"\varphi/\pi",
			ylabel=L"\langle n \rangle",
			#yscale=:log10,
			title=L"\alpha = 2",
			color = [:red :black :green :blue], 
			line = ([2 2 2], [:solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[3, :, r_slider].+var_n_mk[3, :, r_slider] var_n_k[3, :, r_slider] var_n_mk[3, :, r_slider]]),
		),
		leg = :topright,
		#layout = (1, 3),
		#size = (700, 220),
	)
end

# ╔═╡ 9ee646d9-2dc7-481c-82c9-6ad5ccb529ac
data_dict_g2_full["r"][2]

# ╔═╡ bfdb89a9-9e71-41d0-aa87-444c2accb3ae
let 
	ind = 15
	
	a_k_av = (data_dict_g2_full["a_k_av_re"] .+ 1.0im*data_dict_g2_full["a_k_av_im"])
	a_mk_av = (data_dict_g2_full["a_mk_av_re"] .+ 1.0im*data_dict_g2_full["a_mk_av_im"])

	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	var_n_k = variance_n(g2_k, n_k)
	var_n_mk = variance_n(g2_mk, n_mk)

	e_k = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av)
	e_mk = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_mk_av)
	e_full = TI_quantum.e_field_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av, a_mk_av)

	omega = r_list[ind]
	

	fig = plot(
		plot(phi_list/pi, 
			[
				n_k[1, :, ind] .+ n_mk[1, :, ind] .+ 0.5,
				n_k[1, :, ind] .+ 0.5,
				n_mk[1, :, ind] .+ 0.5,
				(2*[alpha_list[1]^2 for j in eachindex(phi_list)] .+ 0.5)/omega,
			].*omega,
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init energy"], 
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\hbar \omega (\langle n \rangle + 1/2)",
			title=L"\alpha = 0.5",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[1, :, ind].+var_n_mk[1, :, ind] var_n_k[1, :, ind] var_n_mk[1, :, ind]]),
		),
		plot(phi_list/pi, 
			[
				n_k[2, :, ind] .+ n_mk[2, :, ind] .+ 0.5,
				n_k[2, :, ind] .+ 0.5,
				n_mk[2, :, ind] .+ 0.5,
				(2*[alpha_list[2]^2 for j in eachindex(phi_list)] .+ 0.5)/omega,
			].*omega,
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init energy"], 
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\hbar \omega (\langle n \rangle + 1/2)",
			title=L"\alpha = 1.25",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[2, :, ind].+var_n_mk[2, :, ind] var_n_k[2, :, ind] var_n_mk[2, :, ind]]),
		),
		plot(phi_list/pi, 
			[
				n_k[3, :, ind] .+ n_mk[3, :, ind] .+ 0.5,
				n_k[3, :, ind] .+ 0.5,
				n_mk[3, :, ind] .+ 0.5,
				(2*[alpha_list[3]^2 for j in eachindex(phi_list)] .+ 0.5)/omega,
			].*omega,
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init energy"], 
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\hbar \omega (\langle n \rangle + 1/2)",
			title=L"\alpha = 2",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[3, :, ind].+var_n_mk[3, :, ind] var_n_k[3, :, ind] var_n_mk[3, :, ind]]),
		),
		#leg = :topright,
		layout = (1, 3),
		size = (700, 220),
		bottom_margin=5.0*Plots.mm,
		left_margin=2.0*Plots.mm,
	)
	NAME_PART = "energy_equal_fields_erasing_"
	savefig(fig, PATH_FIGS*NAME_PART*"r"*string(round(r_list[ind], digits=2))*".pdf")
	fig
end

# ╔═╡ 9fe3e812-0215-4928-a51a-6e0d2543ea72
let 
	ind = 2
	
	a_k_av = (data_dict_g2_full["a_k_av_re"] .+ 1.0im*data_dict_g2_full["a_k_av_im"])
	a_mk_av = (data_dict_g2_full["a_mk_av_re"] .+ 1.0im*data_dict_g2_full["a_mk_av_im"])

	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	alpha_list = data_dict_g2_full["alpha"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	var_n_k = variance_n(g2_k, n_k)
	var_n_mk = variance_n(g2_mk, n_mk)

	e_k = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av)
	e_mk = TI_quantum.e_field_k_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_mk_av)
	e_full = TI_quantum.e_field_average(x_slider, 2*pi, t_slider, 2*pi*r_list[ind], a_k_av, a_mk_av)

	omega = r_list[ind]
	

	fig = plot(
		plot(phi_list/pi, 
			[
				n_k[1, :, ind] .+ n_mk[1, :, ind],
				n_k[1, :, ind],
				n_mk[1, :, ind],
				(2*[alpha_list[1]^2 for j in eachindex(phi_list)]),
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init sum"], 
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 0.5",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[1, :, ind].+var_n_mk[1, :, ind] var_n_k[1, :, ind] var_n_mk[1, :, ind]]),
		),
		plot(phi_list/pi, 
			[
				n_k[2, :, ind] .+ n_mk[2, :, ind],
				n_k[2, :, ind],
				n_mk[2, :, ind],
				(2*[alpha_list[2]^2 for j in eachindex(phi_list)]),
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init sum"], 
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 1.25",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[2, :, ind].+var_n_mk[2, :, ind] var_n_k[2, :, ind] var_n_mk[2, :, ind]]),
		),
		plot(phi_list/pi, 
			[
				n_k[3, :, ind] .+ n_mk[3, :, ind],
				n_k[3, :, ind],
				n_mk[3, :, ind],
				(2*[alpha_list[3]^2 for j in eachindex(phi_list)]),
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k" "init sum"], 
			xlabel=L"\varphi, \; \pi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 2",
			color = [:red :black :black :green], 
			line = ([2 2 2], [:solid :solid :dash :dot]),
			#ribbon = sqrt.([var_n_k[3, :, ind].+var_n_mk[3, :, ind] var_n_k[3, :, ind] var_n_mk[3, :, ind]]),
		),
		#leg = :topright,
		layout = (1, 3),
		size = (700, 220),
		bottom_margin=5.0*Plots.mm,
		left_margin=2.0*Plots.mm,
	)
	NAME_PART = "n_equal_fields_erasing_"
	#savefig(fig, PATH_FIGS*NAME_PART*"r"*string(round(r_list[ind], digits=2))*".pdf")
	fig
end

# ╔═╡ Cell order:
# ╟─e9703909-861e-4c3e-858f-130f7e08d869
# ╠═ed82f1be-e061-11ed-0a27-757649bc8226
# ╠═fef9cfab-5ce3-4bd4-8f2f-fd8ba840e06b
# ╠═8e4f0f92-2284-4408-9ea9-4a2085753643
# ╠═3b52bf51-8e85-4b2f-bfac-a752ee0658ac
# ╠═19b4d18a-0014-486d-9809-d2ac1613393d
# ╠═5080ed27-3963-4a70-80c8-a64870189189
# ╟─bed0ed98-13b5-4e7b-9ffc-f4c52bec3dd0
# ╠═92cac3ba-455f-45f3-a883-c753c053fc96
# ╟─51580c92-f227-4d2a-a378-874ba047294c
# ╠═6cc9021d-2958-4778-acfc-2c189e6454a3
# ╠═68511105-2e94-4bcd-b23a-e2fbc2411f5a
# ╠═7989d94e-a548-4619-adc3-321ff80d0858
# ╟─216f1946-b172-4b0b-b6e2-0f64015f88f9
# ╠═9ee646d9-2dc7-481c-82c9-6ad5ccb529ac
# ╠═bfdb89a9-9e71-41d0-aa87-444c2accb3ae
# ╠═9fe3e812-0215-4928-a51a-6e0d2543ea72
