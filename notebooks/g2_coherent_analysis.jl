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

	NMAX = 50
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

# ╔═╡ ec83206d-1088-4ac7-b11f-a37357a903cf
# ╠═╡ disabled = true
#=╠═╡
plot(

		# before switching	
		plot(x_list, 
			[
				[real(e_field_average(i, k, t_slider, omega, alpha, beta*exp(im*phi_list[1]))) for i in x_list], 
				[real(e_field_k_average(i, k, t_slider, omega, alpha)) for i in x_list],
				[real(e_field_k_average(i, -k, t_slider, omega, beta*exp(im*phi_list[1]))) for i in x_list],
			],
			ylims=(-20,20), 
			label=["full" "k" "-k"], 
			title=L"\Delta \varphi = 0",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),

		plot(x_list, 
			[
				[real(e_field_average(i, k, t_slider, omega, alpha, beta*exp(im*phi_list[2]))) for i in x_list], 
				[real(e_field_k_average(i, k, t_slider, omega, alpha)) for i in x_list],
				[real(e_field_k_average(i, -k, t_slider, omega, beta*exp(im*phi_list[2]))) for i in x_list],
			],
			ylims=(-20,20), 
			label=["full" "k" "-k"], 
			title=L"\Delta \varphi= \pi /2",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),

		plot(x_list, 
			[
				[real(e_field_average(i, k, t_slider, omega, alpha, beta*exp(im*phi_list[3]))) for i in x_list], 
				[real(e_field_k_average(i, k, t_slider, omega, alpha)) for i in x_list],
				[real(e_field_k_average(i, -k, t_slider, omega, beta*exp(im*phi_list[3]))) for i in x_list],
			],
			ylims=(-20,20), 
			label=["full" "k" "-k"], 
			title=L"\Delta \varphi= \pi",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),

		
		# after switching
		plot(x_list, 
			[
				[real(e_field_average(i, 2*pi, t_slider, 2*pi/r_i_slider, a_k_av[1, r_i_slider], a_mk_av[1, r_i_slider])) for i in x_list],
				[real(e_field_k_average(i, 2*pi, t_slider, 2*pi/r_i_slider, a_k_av[1, r_i_slider])) for i in x_list],
				[real(e_field_k_average(i, -2*pi, t_slider, 2*pi/r_i_slider, a_mk_av[1, r_i_slider])) for i in x_list]
			],
			ylims=(-20,20), 
			label=["full" "k" "-k"], 
			#title="after switching",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(x_list, 
			[
				[real(e_field_average(i, 2*pi, t_slider, 2*pi/r_i_slider, a_k_av[2, r_i_slider], a_mk_av[2, r_i_slider])) for i in x_list],
				[real(e_field_k_average(i, 2*pi, t_slider, 2*pi/r_i_slider, a_k_av[2, r_i_slider])) for i in x_list],
				[real(e_field_k_average(i, -2*pi, t_slider, 2*pi/r_i_slider, a_mk_av[2, r_i_slider])) for i in x_list]
			],
			ylims=(-20,20), 
			label=["full" "k" "-k"], 
			title="after switching",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		plot(x_list, 
			[
				[real(e_field_average(i, 2*pi, t_slider, 2*pi/r_i_slider, a_k_av[3, r_i_slider], a_mk_av[3, r_i_slider])) for i in x_list],
				[real(e_field_k_average(i, 2*pi, t_slider, 2*pi/r_i_slider, a_k_av[3, r_i_slider])) for i in x_list],
				[real(e_field_k_average(i, -2*pi, t_slider, 2*pi/r_i_slider, a_mk_av[3, r_i_slider])) for i in x_list]
			],
			ylims=(-20,20), 
			label=["full" "k" "-k"], 
			#title="after switching",
			color = [:red :black :black], 
			line = ([2 1 1], [:solid :dash :dot]),
		),
		leg = :topright,
		size = (800, 600),
	)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─e9703909-861e-4c3e-858f-130f7e08d869
# ╠═ed82f1be-e061-11ed-0a27-757649bc8226
# ╠═fef9cfab-5ce3-4bd4-8f2f-fd8ba840e06b
# ╠═025d67a9-ce1e-47ab-ae28-b953573440c1
# ╠═21d91c6d-df88-4b6e-b598-e2192fa6114e
# ╠═74e1988c-705f-4881-8606-6c7e3712174f
# ╠═a0e9722f-fa78-494f-b380-d8cc2d2c3ffc
# ╠═ec83206d-1088-4ac7-b11f-a37357a903cf
