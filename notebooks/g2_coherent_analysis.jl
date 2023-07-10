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

# ╔═╡ 025d67a9-ce1e-47ab-ae28-b953573440c1
begin
	# Load data
	alpha_list = [0.5, 2.0]
	data_dict_g2 = []
	for alpha in alpha_list
		NAME_PART = "alpha"*string(alpha)*"_beta_r_dependent_more_points"*".h5"
		push!(data_dict_g2, load(PATH_DATA*"g2_erasing_"*NAME_PART))
	end
end

# ╔═╡ 21d91c6d-df88-4b6e-b598-e2192fa6114e
data_dict_g2[1]

# ╔═╡ 74e1988c-705f-4881-8606-6c7e3712174f
@bind phi_slider PlutoUI.Slider(1:3, default=1)

# ╔═╡ a0e9722f-fa78-494f-b380-d8cc2d2c3ffc
let
	phi_list = data_dict_g2[1]["phi"]
	r_list = data_dict_g2[1]["r"]
	beta_list = [data_dict_g2[i]["beta"]  for i in eachindex(alpha_list)]
	g2_k = [data_dict_g2[i]["g2_k"]  for i in eachindex(alpha_list)]
	g2_mk = [data_dict_g2[i]["g2_mk"]  for i in eachindex(alpha_list)]
	g2_cross = [data_dict_g2[i]["g2_cross"]  for i in eachindex(alpha_list)]

	plot(
		plot(r_list, 
			[
				g2_cross[1][phi_slider, :],
				g2_k[1][phi_slider, :],
				g2_mk[1][phi_slider, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 0.5",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				g2_cross[2][phi_slider, :],
				g2_k[2][phi_slider, :],
				g2_mk[2][phi_slider, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 2",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		leg = :topright,
		#size = (800, 600),
	)
end

# ╔═╡ 3b52bf51-8e85-4b2f-bfac-a752ee0658ac
begin
	# Load data
	NAME_PART_full = "alpha"*"_beta_r_dependent_more_points"*".h5"
	data_dict_g2_full = load(PATH_DATA*"g2_field_erasing_"*NAME_PART_full)
end

# ╔═╡ 2eba5b39-0f79-4555-baa7-383acd4498bb
@bind phi_slider_full PlutoUI.Slider(1:9, default=1)

# ╔═╡ bed0ed98-13b5-4e7b-9ffc-f4c52bec3dd0
let
	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	plot(
		plot(r_list, 
			[
				g2_cross[1, phi_slider_full, :],
				g2_k[1, phi_slider_full, :],
				g2_mk[1, phi_slider_full, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 0.5",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				g2_cross[2, phi_slider_full, :],
				g2_k[2, phi_slider_full, :],
				g2_mk[2, phi_slider_full, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 1.25",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				g2_cross[3, phi_slider_full, :],
				g2_k[3, phi_slider_full, :],
				g2_mk[3, phi_slider_full, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 2",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[1, phi_slider_full, :] .+ n_mk[1, phi_slider_full, :],
				n_k[1, phi_slider_full, :],
				n_mk[1, phi_slider_full, :]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			#title=L"\alpha = 0.5",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[2, phi_slider_full, :] .+ n_mk[2, phi_slider_full, :],
				n_k[2, phi_slider_full, :],
				n_mk[2, phi_slider_full, :]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			#title=L"\alpha = 1.25",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[3, phi_slider_full, :] .+ n_mk[3, phi_slider_full, :],
				n_k[3, phi_slider_full, :],
				n_mk[3, phi_slider_full, :]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			#title=L"\alpha = 2",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[1, phi_slider_full, :], n_k[1, phi_slider_full, :]),
				variance_n(g2_mk[1, phi_slider_full, :], n_mk[1, phi_slider_full, :])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 0.5",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[2, phi_slider_full, :], n_k[2, phi_slider_full, :]),
				variance_n(g2_mk[2, phi_slider_full, :], n_mk[2, phi_slider_full, :])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 1.25",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[3, phi_slider_full, :], n_k[3, phi_slider_full, :]),
				variance_n(g2_mk[3, phi_slider_full, :], n_mk[3, phi_slider_full, :])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 2",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		leg = :topright,
		#size = (800, 600),
	)
end

# ╔═╡ 92cac3ba-455f-45f3-a883-c753c053fc96
let

	ind = 1 #  in the range 1:9 -- 0:2π
	
	phi_list = data_dict_g2_full["phi"]
	r_list = data_dict_g2_full["r"]
	beta_list = data_dict_g2_full["beta"]
	
	g2_k = data_dict_g2_full["g2_k"]
	g2_mk = data_dict_g2_full["g2_mk"]
	g2_cross = data_dict_g2_full["g2_cross"]

	n_k = data_dict_g2_full["n_k"]
	n_mk = data_dict_g2_full["n_mk"]

	fig = plot(
		plot(r_list, 
			[
				g2_cross[1, ind, :],
				g2_k[1, ind, :],
				g2_mk[1, ind, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 0.5",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				g2_cross[2, ind, :],
				g2_k[2, ind, :],
				g2_mk[2, ind, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 1.25",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				g2_cross[3, ind, :],
				g2_k[3, ind, :],
				g2_mk[3, ind, :]
			],
			#ylims=(0,2), 
			label=["cross" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"g^{(2)}(0)",
			title=L"\alpha = 2",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[1, ind, :] .+ n_mk[1, ind, :],
				n_k[1, ind, :],
				n_mk[1, ind, :]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			#title=L"\alpha = 0.5",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[2, ind, :] .+ n_mk[2, ind, :],
				n_k[2, ind, :],
				n_mk[2, ind, :]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			#title=L"\alpha = 1.25",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				n_k[3, ind, :] .+ n_mk[3, ind, :],
				n_k[3, ind, :],
				n_mk[3, ind, :]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\langle n \rangle",
			#title=L"\alpha = 2",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[1, ind, :], n_k[1, ind, :]),
				variance_n(g2_mk[1, ind, :], n_mk[1, ind, :])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 0.5",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[2, ind, :], n_k[2, ind, :]),
				variance_n(g2_mk[2, ind, :], n_mk[2, ind, :])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 1.25",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		plot(r_list, 
			[
				variance_n(g2_k[3, ind, :], n_k[3, ind, :]),
				variance_n(g2_mk[3, ind, :], n_mk[3, ind, :])
			],
			#ylims=(0,2), 
			label=["k" "-k"], 
			xlabel=L"r = n_0/n_1",
			ylabel=L"\mathrm{Var}(n)",
			#title=L"\alpha = 2",
			color = [:black :black], 
			line = ([1 1], [:dash :dot]),
		),
		leg = :topright,
		#size = (800, 600),
	)
	NAME_PART = "g2_erasing_"
	#savefig(fig, PATH_FIGS*NAME_PART*"phi"*string(phi_list[ind])*".pdf")
	fig
end

# ╔═╡ 51580c92-f227-4d2a-a378-874ba047294c
md"""
 $r$:  $(@bind r_slider PlutoUI.Slider(1:NMAX, default=NMAX))

 $t$:  $(@bind t_slider PlutoUI.Slider(0:400, default=0))

 $x$:  $(@bind x_slider PlutoUI.Slider(0:400, default=0))

 $\alpha$:  $(@bind vphi_slider PlutoUI.Slider(1:3, default=1))
"""

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
		plot(phi_list, 
			[
				n_k[1, :, r_slider] .+ n_mk[1, :, r_slider],
				n_k[1, :, r_slider],
				n_mk[1, :, r_slider]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"\varphi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 0.5",
			color = [:red :black :green], 
			line = ([2 2 2], [:solid :dash :dot]),
			ribbon = sqrt.([var_n_k[1, :, r_slider].+var_n_mk[1, :, r_slider] var_n_k[1, :, r_slider] var_n_mk[1, :, r_slider]]),
		),
		plot(phi_list, 
			[
				n_k[2, :, r_slider] .+ n_mk[2, :, r_slider],
				n_k[2, :, r_slider],
				n_mk[2, :, r_slider]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"\varphi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 1.25",
			color = [:red :black :green], 
			line = ([2 2 2], [:solid :dash :dot]),
			ribbon = sqrt.([var_n_k[2, :, r_slider].+var_n_mk[2, :, r_slider] var_n_k[2, :, r_slider] var_n_mk[2, :, r_slider]]),
		),
		plot(phi_list, 
			[
				n_k[3, :, r_slider] .+ n_mk[3, :, r_slider],
				n_k[3, :, r_slider],
				n_mk[3, :, r_slider]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"\varphi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 2",
			color = [:red :black :green], 
			line = ([2 2 2], [:solid :dash :dot]),
			ribbon = sqrt.([var_n_k[3, :, r_slider].+var_n_mk[3, :, r_slider] var_n_k[3, :, r_slider] var_n_mk[3, :, r_slider]]),
		),
		leg = :topright,
		#layout = (1, 3),
		#size = (700, 220),
	)
end

# ╔═╡ bfdb89a9-9e71-41d0-aa87-444c2accb3ae
let 
	ind = 1
	
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

	fig = plot(
		plot(phi_list, 
			[
				n_k[1, :, ind] .+ n_mk[1, :, ind],
				n_k[1, :, ind],
				n_mk[1, :, ind]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"\varphi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 0.5",
			color = [:red :black :green], 
			line = ([2 2 2], [:solid :dash :dot]),
			ribbon = sqrt.([var_n_k[1, :, ind].+var_n_mk[1, :, ind] var_n_k[1, :, ind] var_n_mk[1, :, ind]]),
		),
		plot(phi_list, 
			[
				n_k[2, :, ind] .+ n_mk[2, :, ind],
				n_k[2, :, ind],
				n_mk[2, :, ind]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"\varphi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 1.25",
			color = [:red :black :green], 
			line = ([2 2 2], [:solid :dash :dot]),
			ribbon = sqrt.([var_n_k[2, :, ind].+var_n_mk[2, :, ind] var_n_k[2, :, ind] var_n_mk[2, :, ind]]),
		),
		plot(phi_list, 
			[
				n_k[3, :, ind] .+ n_mk[3, :, ind],
				n_k[3, :, ind],
				n_mk[3, :, ind]
			],
			#ylims=(0,2), 
			label=["sum" "k" "-k"], 
			xlabel=L"\varphi",
			ylabel=L"\langle n \rangle",
			title=L"\alpha = 2",
			color = [:red :black :green], 
			line = ([2 2 2], [:solid :dash :dot]),
			ribbon = sqrt.([var_n_k[3, :, ind].+var_n_mk[3, :, ind] var_n_k[3, :, ind] var_n_mk[3, :, ind]]),
		),
		leg = :topright,
		layout = (1, 3),
		size = (700, 220),
		bottom_margin=5.0*Plots.mm,
	)
	NAME_PART = "n_varn_erasing_"
	#savefig(fig, PATH_FIGS*NAME_PART*"r"*string(r_list[ind])*".pdf")
	fig
end

# ╔═╡ Cell order:
# ╟─e9703909-861e-4c3e-858f-130f7e08d869
# ╠═ed82f1be-e061-11ed-0a27-757649bc8226
# ╠═fef9cfab-5ce3-4bd4-8f2f-fd8ba840e06b
# ╠═8e4f0f92-2284-4408-9ea9-4a2085753643
# ╠═025d67a9-ce1e-47ab-ae28-b953573440c1
# ╠═21d91c6d-df88-4b6e-b598-e2192fa6114e
# ╟─74e1988c-705f-4881-8606-6c7e3712174f
# ╟─a0e9722f-fa78-494f-b380-d8cc2d2c3ffc
# ╟─3b52bf51-8e85-4b2f-bfac-a752ee0658ac
# ╠═2eba5b39-0f79-4555-baa7-383acd4498bb
# ╠═bed0ed98-13b5-4e7b-9ffc-f4c52bec3dd0
# ╠═92cac3ba-455f-45f3-a883-c753c053fc96
# ╟─51580c92-f227-4d2a-a378-874ba047294c
# ╟─216f1946-b172-4b0b-b6e2-0f64015f88f9
# ╠═bfdb89a9-9e71-41d0-aa87-444c2accb3ae
