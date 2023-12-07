### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ 9a8f509a-e8ce-4324-9c78-7fe13b4d3af1
begin
	using Plots
	using BenchmarkTools
	using PlutoUI
	using HDF5, FileIO, Printf, LaTeXStrings
	using LinearAlgebra

	# PATH_FIGS, PATH_DATA = path()
	PATH_FIGS, PATH_DATA = ("../../Figs/", "../../Data/")
end

# ╔═╡ 08c2d5fe-e42c-416a-8eb5-90ed956056d3
begin
	using PGFPlotsX
end

# ╔═╡ 12f49f28-7d82-11ee-18e0-8f83344ebab3
md"
# Data analysis: Time Interfaces
"

# ╔═╡ c5bdf327-a439-4b3e-95c0-de3f8d45c970
"""
    path()

Returns the paths for data and figures, depending on the home directory of the user. If the home directory is "C:\\Users\\nnefedkin", it returns the paths for Windows, if it's "/home/nikita", it returns the paths for Linux, and if it's "/Users/jimi", it returns the paths for MacOS.

# Returns:
- PATH_FIGS: The path for the figures.
- PATH_DATA: The path for the data.
"""
function path()
    home = homedir()
    if home == "C:\\Users\\nnefedkin"
        PATH_FIGS = "D:/nnefedkin/Google_Drive/Work/In process/Projects/Time_interfaces/Figs/"
        PATH_DATA = "D:/nnefedkin/Google_Drive/Work/In process/Projects/Time_interfaces/Data/"
    elseif home == "/home/nikita"
        PATH_FIGS = "/home/nikita/Documents/figs/time_interfaces/"
        PATH_DATA = "/home/nikita/Documents/data/time_interfaces/"
    elseif home == "/Users/jimi"
        PATH_FIGS = "/Users/jimi/Google Drive/Work/In process/Projects/Time_interfaces/Figs/"
        PATH_DATA = "/Users/jimi/Google Drive/Work/In process/Projects/Time_interfaces/Data/"
    end
    return PATH_FIGS, PATH_DATA
end

# ╔═╡ 57b440e0-3324-4a03-b2a1-e1b8cdfb986a
begin
	function n_photons_vacuum(u, v)
		return 2*(1 - abs(v/u)^2)*(u*v)^2 / (u^2 - v^2)^2
	end
	function u(r)
		return 0.5 * (1 + r) / sqrt(r)
	end
	function v(r)
		return 0.5 * (1 - r) / sqrt(r)
	end
	function ent_entropy_vac(r)
		return (- (1-r)^2 * log((1-r^2/(1+r)^2)) / 4 / r - log(4*r/(1+r)^2))
	end
end

# ╔═╡ c5f0ba7a-8b39-4f26-8579-232f95a05279
"""
    load_dict(data_type::String; alpha::Real = 2.0, beta::Real = 2.0, numsum::Int = 100, mode::String = "equal", r_range::String = "")

Load a dictionary from an HDF5 file based on the specified data type and parameters.

# Arguments

- `data_type::String`: The type of data to load. Supported values are "coh_photon_distrib," "g2," "operators," "entropy."

# Keyword Arguments

- `alpha::Real = 2.0`: Amplitude of k-mode (default: 2.0).
- `beta::Real = 2.0`: Amplitude of -k-mode (default: 2.0).
- `numsum::Int = 100`: Number of terms in the sum (default: 100).
- `mode::String = "equal"`: if alpha/beta is "equal" or in the "erasing" mode (default: "equal").
- `r_range::String = ""`: "rgg1_" if r > 1, or "" if r < 1 (default: "").

# Returns

- A dictionary loaded from the corresponding HDF5 file.

The function constructs the file name based on the `data_type` and the provided parameters, and then loads the dictionary from the file located at `PATH_DATA * FILE_NAME`.

# Examples

```julia
dict = load_dict("coh_photon_distrib"; alpha=3.0, beta=3.0, numsum=100)
```
"""
function load_dict(data_type::String; alpha::Real = 2.0, beta::Real = 2.0,
				   numsum::Int = 100, mode::String = "equal", r_range::String = "")
	if data_type == "coh_photon_distrib"
		FILE_NAME = data_type*"_alpha"*string(alpha)*"_beta"*string(beta)*
					"_numsum"*string(numsum)*".h5"
	elseif data_type == "coh_photon_distrib_more_points"
		FILE_NAME = "coh_photon_distrib_alpha"*string(alpha)*"_beta"*string(beta)*
					"_numsum"*string(numsum)*"_more_points.h5"
	elseif data_type == "g2"
		if mode == "erasing"
			FILE_NAME = data_type*"_field_"*mode*"_alpha_beta_r_dependent"*
						"_more_points"*".h5"
		elseif mode == "equal"
			FILE_NAME = data_type*"_field_"*mode*"_phi_"*r_range*
						"alpha_beta_r_dependent"*"_more_points"*".h5"
		end
	elseif data_type == "operators"
		if mode == "equal"
			FILE_NAME = data_type*"_average_"*"alpha"*string(Int(alpha))*
						"_beta"*string(Int(beta))*".h5"
		elseif mode == "erasing"
			FILE_NAME = data_type*"_average_"*mode*"_alpha"*string(Int(alpha))*
						"_beta_r_dependent"*".h5"
		end
	elseif data_type == "entropy"
		FILE_NAME = data_type*"_ent_nterms"*string(numsum)*"_"*mode*
					"_alpha_beta"*".h5"
	end
	return load(PATH_DATA*FILE_NAME)
end


# ╔═╡ 766af561-491d-40d9-92db-80fd2956bfa4
dict = load_dict("g2"; alpha = 2.0, mode = "equal", r_range="rgg1_")

# ╔═╡ a7eaa2b8-0b68-4987-8d3f-334595cf1891
# Sets of avaliable parameters
begin
	data_types = ["coh_photon_distrib", "g2", "operators"]
	numsums = [100, 200]
	alphas = [0.5, 2.0, 3.0, 5.0, 10.0]
	betas = [3.0, 5.0, 10.0]
	modes = ["erasing", "equal"]
	r_range = ["", "rgg1_"]
end

# ╔═╡ 50ae134f-133a-4f65-8cf4-ff2256300d54
# Load data
begin
	data_photon_distrib = Dict([("alpha_beta="*string(i), 
								 load_dict("coh_photon_distrib";
								 alpha = i, beta = i,
								 numsum = (i == 10.0) ? 200 : 100)) 
								for i in betas])
	data_photon_distrib_mp = Dict([("alpha_beta="*string(5.0),
								 load_dict("coh_photon_distrib_more_points";
								 alpha = 5.0, beta = 5.0,
								 numsum = 100))])
	data_g2 = Dict([(i == "erasing") ? 
					(i, load_dict("g2"; mode = i)) : 
					(i, Dict([("rll1", load_dict("g2"; mode = i)), 
							  ("rgg1", load_dict("g2"; mode = i, r_range="rgg1_"))]))
					for i in modes])
	data_operators = Dict([(i, 
							load_dict("operators"; mode = i, alpha=2.0, beta=2.0))
						  for i in modes])
	data_entropy = Dict([("equal", 
						  load_dict("entropy"; mode = "equal", numsum = 30))])
end

# ╔═╡ 52cef578-5405-49a1-853d-58297fd15e25
data_g2

# ╔═╡ 63078392-10a0-458d-b80b-d6fc3e34dfcf
let
	x = range(0.0, pi, 9)
	y = data_entropy["equal"]["r"]
	z = data_entropy["equal"]["entropy_ent_vn"]

	p = plot(zeros(length(y)) .+ x[1], y, z[1,1,:])
	for i in 2:length(x)
		plot!(zeros(length(y)) .+ x[i], y, z[1,i,:])
	end
	p
end

# ╔═╡ 308bc872-f072-4399-815b-3a2fb562a2d4
maximum(data_entropy["equal"]["entropy_ent_vn"])

# ╔═╡ 3147dee6-4685-42a8-89c5-79028e8cc3a6
md"
## Energy and number of photons

**Variance**: $(g^{(2)}(0) - 1) \langle a^\dagger a \rangle^2 + \langle a^\dagger a \rangle$

Idea: competition between vacuum photons and initial photons at ti
"

# ╔═╡ 0f9f6870-b5a7-4461-b58f-071889d7bc33
data_g2["equal"]["rgg1"]

# ╔═╡ 2a25376b-ee9f-4172-80ae-fba2c744f70f
data_g2["equal"]["rll1"]["alpha"] .* data_g2["equal"]["rll1"]["n_k"][:,:,20]

# ╔═╡ 15b0c4fc-7765-4036-b42a-982ec8886cc9
begin
	slider_r = @bind r_slider PlutoUI.Slider(eachindex(data_g2["equal"]["rll1"]["r"]), default=1);
	slider_phi = @bind phi_slider PlutoUI.Slider(eachindex(data_g2["equal"]["rll1"]["phi"]), default=1);
	slider_alpha = @bind alpha_slider PlutoUI.Slider(eachindex(data_g2["equal"]["rll1"]["alpha"]), default=1);
end;

# ╔═╡ 79439794-d465-4ec8-ab81-27ef04110f51
md"
$slider_r r
"

# ╔═╡ 8320039b-776b-493d-a2d0-4ee56dc0d7a5
print("r = ", round(data_g2["equal"]["rll1"]["r"][r_slider], digits=3), ", ϕ = ", data_g2["equal"]["rll1"]["phi"][phi_slider] / π, "π")

# ╔═╡ 1c1ce53f-03c8-4ae3-b40f-6060b558e953
# equal 
let
	gr()

	alphas = data_g2["equal"]["rll1"]["alpha"]
	rs = data_g2["equal"]["rll1"]["r"]
	r_entropy = range(0.01, 1, 100)
	
	omega_0 = 1.0
	omega_1 = omega_0 * rs[r_slider]
	init_energy = omega_0 .* (2*alphas.^2 .+ 0.5)
	init_photons = 2*alphas.^2
	
	
	x = data_g2["equal"]["rll1"]["phi"]
	y_1 = (data_g2["equal"]["rll1"]["n_k"][:,:,r_slider] .+ 
			data_g2["equal"]["rll1"]["n_mk"][:,:,r_slider]) # number of photons
	y_2 = omega_1 .* ((data_g2["equal"]["rll1"]["n_k"][:,:,r_slider] .+ 
			data_g2["equal"]["rll1"]["n_mk"][:,:,r_slider]) .+ 0.5) ./ init_energy # energy
	y_3 = ((data_g2["equal"]["rll1"]["g2_k"][:,:,r_slider] .- 1.0) .* 
			data_g2["equal"]["rll1"]["n_k"][:,:,r_slider].^2 .+ 
			data_g2["equal"]["rll1"]["n_k"][:,:,r_slider])./(y_1./2) # variance in 1 mode
	y_4 = data_entropy["equal"]["entropy_ent_vn"][1,1,:]
	y_4_vac = ent_entropy_vac.(r_entropy)

	l = @layout([grid(2,1){0.6w} grid(3,1)])

	kw = (
		# xlabel=L"\Delta \varphi",
		xticks=([0, π/2, π], ["0", "π/2", "π"]),
		label=[L"|\alpha| = 0.5" L"|\alpha| = 2"],
	)
	
	p_1 = plot(x, [y_2'[:,1] y_2'[:,3]], 
		lw=2, m=4,
		ls=[:solid :dash],
		ylims=(-0.05, 1.05),
		# yscale=:log10,
		ylabel=L"$E / E_0$";
		kw...
	)
	hline!([1], c=:black, ls=:dot, label=:none)
	p_2 = plot(x, [y_1'[:,1] y_1'[:,3]], 
		lw=2, m=4,
		ls=[:solid :dash],
		ylims=(-0.1, maximum(y_1'[:,3]) + 2),
		# yscale=:log10,
		ylabel=L"$\langle n \rangle$";
		kw...
	)
	p_3 = plot(x, [y_3'[:,1] y_3'[:,3]],
		lw=2, m=4,
		ls=[:solid :dash],
		legend=:bottomleft,
		# yscale=:log10,
		xlabel=L"\Delta \varphi",
		ylabel=L"$\mathrm{Var}(n) / \langle n \rangle$";
		kw...
	)
	p_4 = plot(r_entropy, y_4_vac,
		lw=2, 
		xlabel=L"r = n_0 / n_1",
		ylabel=L"S_\mathrm{ent}",
		legend=false,
	)
	scatter!([rs[r_slider]], [ent_entropy_vac(rs[r_slider])], m=6, c=:red)
	p_5 = plot([0], [0])

	p = plot(p_5, p_4, p_1, p_2, p_3, layout=l, size=(500,500))
	# savefig(p, PATH_FIGS*"Fig_1_with_empty_v1.1.pdf")
end

# ╔═╡ 5c157c15-ad97-4a1a-bec5-81975118a8ab
md"
$slider_phi ϕ

$slider_alpha α
"

# ╔═╡ 9212af44-2ce8-4459-bab2-b2aabfc2b1a4
print("ϕ = ", data_g2["equal"]["rll1"]["phi"][phi_slider] / π, "π", ", α = ", data_g2["equal"]["rll1"]["alpha"][alpha_slider])

# ╔═╡ 9c70ceef-d560-4d2b-b791-27dae6eada90
# equal 
let
	gr()

	alphas = data_g2["equal"]["rll1"]["alpha"]
	rs = data_g2["equal"]["rll1"]["r"][1:end]
	
	omega_0 = 1.0
	omega_1 = omega_0 * rs[r_slider]
	init_energy = omega_0 .* (2*alphas.^2 .+ 0.5)
	init_photons = 2*alphas.^2

	n_k = data_g2["equal"]["rll1"]["n_k"][alpha_slider,phi_slider,1:end]
	n_mk = data_g2["equal"]["rll1"]["n_mk"][alpha_slider,phi_slider,1:end]
	g2_k = data_g2["equal"]["rll1"]["g2_k"][alpha_slider,phi_slider,1:end]
	g2_cross = data_g2["equal"]["rll1"]["g2_cross"][alpha_slider,phi_slider,1:end]
	var_nk = (g2_k .- 1.0) .* n_k.^2 .+ n_k
 	
	
	x = rs
	y_1 = (n_k .+ n_mk) # number of photons
	y_1_class = 2*alphas[alpha_slider]^2 .* x
	y_1_vac = n_photons_vacuum.(u.(x), v.(x))
	y_2_0 = g2_k
	y_2_1 = g2_cross
	y_3 = var_nk # variance in 1 mode

	l = @layout([a; b; c])

	kw = (
		# xlabel=L"r = n_0 / n_1",
		leftmargin=3Plots.mm,
	)

	p_1 = plot(x, [y_1 y_1_class y_1_vac], 
		lw=2,
		ls=[:solid :dash :solid],
		# yscale=:log10,
		ylabel=L"$\langle n \rangle$",
		label=["quantum" "classical" "vacuum" "q+c"];
		kw...
	)
	p_2 = plot(x, [y_2_0, y_2_1], 
		lw=2,
		ls=[:solid :dash],
		# yscale=:log10,
		ylabel=L"g^{(2)}(0)",
		label=[L"k, -k" "cross"];
		kw...
	)
	p_3 = plot(x, y_3 ./ n_k,
		lw=2,
		# yscale=:log10,
		xlabel=L"r = \sqrt{\varepsilon_0 / \varepsilon_1}",
		ylabel=L"\mathrm{Var}(n) / \langle n \rangle",
		label=L"k, -k";
		kw...
	)

	p = plot(p_1, p_2, p_3, layout=l, size=(400,700))
	# savefig(p, PATH_FIGS*"Fig_2_phipi_col.pdf")
end

# ╔═╡ 44500e2a-4eb7-401a-bae0-b2008b7751f2
md"""
## Distributions
"""

# ╔═╡ 1e8a8993-3e0a-4b21-b2b2-8fa85d8994a3
data_photon_distrib

# ╔═╡ e30ab803-024c-47b4-9936-d2dfa82e4dae
begin
	dict_idx_global = "alpha_beta=5.0"
	slider_r_ph = @bind r_slider_ph PlutoUI.Slider(eachindex(data_photon_distrib[dict_idx_global]["r"]), default=1);
	slider_slice = @bind slice_slider PlutoUI.Slider(1:(2*length(data_photon_distrib[dict_idx_global]["n_1"])-3), default=1);
end;

# ╔═╡ 0b0f4ffc-f83f-47e0-8bb6-c99ccdc256ab
md"
$slider_r_ph r (distribution)

$slider_slice slice
"

# ╔═╡ 47c3e694-9e44-4adb-9170-1b18ca10a2f4
print("r = ", round(data_photon_distrib["alpha_beta=10.0"]["r"][r_slider_ph], digits=3), " slice = ", slice_slider)

# ╔═╡ e065559e-d2ff-4734-a8e0-11d69c58f01a
# 2D distributions
let
	gr()
	# pgfplotsx()
	
	dict_idx = dict_idx_global
	rs = data_photon_distrib[dict_idx]["r"]
	r= data_photon_distrib[dict_idx]["r"][r_slider_ph]
	n_av_vac = (1 - r)^2 / (4 * r) 

	
	x = data_photon_distrib[dict_idx]["n_1"]
	y = data_photon_distrib[dict_idx]["n_2"]
	z = data_photon_distrib[dict_idx]["p_ab"][:,:,r_slider_ph]
	z_d = diag(z)
	p_vac = [4*r*(1-r)^(2*m) / (1+r)^(2*(m+1)) for m in x]
	p_vac_m = diagm(p_vac)

	x_l = (slice_slider ÷ length(x) < 1) ? [x[end]-slice_slider, x[end]] : [x[1], 2*x[end] - slice_slider]
	y_l = (slice_slider ÷ length(x) < 1) ? [y[end], y[end]-slice_slider] : [2*y[end]-slice_slider, y[1]]
	z_b = ((slice_slider ÷ length(x) < 1) ?
			[z[length(x)-i, length(y)-slice_slider+i] for i in 0:slice_slider] :
			[z[length(x)-(slice_slider%length(x)+1)-i, 1 + i] for i in 0:(length(x) - (slice_slider%length(x)+2))])

	kw = (
		xlabel = L"N_k",
		ylabel = L"N_{-k}",
	)

	l = @layout([
		a{0.6w} grid(2,1)
	])
	
	p_1 = heatmap(x, y, (z), 
		aspect_ratio=1;
		kw...
	)
	plot!(x_l, y_l, 
		lw=2, c=:white,	
		label=:none,
	)
	scatter!([n_av_vac], [n_av_vac],
		m=5, c=:white,
		label=L"\langle n_{k,-k} \rangle_\mathrm{vac}",
	)
	# heatmap(x, y, p_vac_m, 
	# 	aspect_ratio=1, alpha=0.4;
	# 	kw...
	# )
	p_2 = plot(x, z_d, seriestype=:bar, lw=0, label=:none)
	p_3 = plot(z_b,
		seriestype=:bar,
		lw=0,
		xlabel = L"| N_{k}, N_{-k}\rangle",
		ylabel = L"P(N_k, N_{-k})",
		xrotation = 60,
		label=:none,
		# xticks=(ticks_num, ticks_let),
		# title=L"r = %$labl";
		# kw...
	)
	p = plot(p_1, p_2, p_3, layout=l)
	# print([(length(x)-(slice_slider%length(x)+1)-i, 1 + i) for i in 0:(length(x) - (slice_slider%length(x)+2))])
end

# ╔═╡ 3da71126-bc76-43e9-9004-37a8f37d74d5
let
	gr()
	# pgfplotsx()
	
	dict_idx = "alpha_beta=3.0"
	rs = data_photon_distrib[dict_idx]["r"]

	
	x = data_photon_distrib[dict_idx]["n_1"]
	y = data_photon_distrib[dict_idx]["n_2"]

	kw = (
		xlabel = L"N_k",
		ylabel = L"N_{-k}",
	)
	
 	anim = @animate for i ∈ length(rs):-1:1
		r = data_photon_distrib[dict_idx]["r"][i]
		n_av_vac = (1 - r)^2 / (4 * r) 
		z = data_photon_distrib[dict_idx]["p_ab"][:,:,i]
		heatmap(x, y, z, 
			title = latexstring("r = $(round(r, digits=2))"),
			aspect_ratio=1;
			kw...
		)
		scatter!([n_av_vac], [n_av_vac],
			m=5, c=:white,
			label=L"\langle n_{k,-k} \rangle_\mathrm{vac}",
		)
	end
	gif(anim, PATH_FIGS*"anim_dist_"*string(dict_idx)*"_fps15.gif", fps = 10)
end

# ╔═╡ ef6bb4ad-5ef4-47ef-ba88-f07fae9c3a1a
# 2D distributions 
let
	# gr()
	pgfplotsx()

	idx = [30, 16, 6, 2]
	
	dict_idx = "alpha_beta=5.0"
	rs = data_photon_distrib[dict_idx]["r"]
	r = [data_photon_distrib[dict_idx]["r"][i] for i in idx]
	n_av_vac = (1 .- r).^2 ./ (4 .* r) 

	
	x = data_photon_distrib[dict_idx]["n_1"]
	y = data_photon_distrib[dict_idx]["n_2"]
	z = [data_photon_distrib[dict_idx]["p_ab"][:,:,i] for i in idx]

	kw = (
		xlabel = L"N_k",
		ylabel = L"N_{-k}",
		# framestyle=:box,
		colorbar=:right,
		colorbar_formatter=:scientific,
		titlefontsize=12,
		guidefontsize=12,
		tickfontsize=10,
		colorbar_tickfontsize=10,
		legend=:bottomright,
		colorbar_title = L"P(N_k, N_{-k})",
		colorbar_titlefontsize=10,
		rightmargin=15Plots.mm,
		leftmargin=3Plots.mm,
		bottommargin=8Plots.mm,
	)
	l = @layout([a b c d])

	p_i = []
	for i in eachindex(idx)
		push!(p_i, heatmap(x, y, z[i], 
			title = latexstring("r = $(round(r[i], digits=2))"),
			aspect_ratio=:equal;
			kw...
		))
		scatter!([n_av_vac[i]], [n_av_vac[i]],
			m=4, c=:white,
			label=L"\langle n_{k,-k} \rangle_\mathrm{vac}",
		)
	end
	p = plot(p_i[1], p_i[2], p_i[3], p_i[4], layout=l, size=(1500,300))
	savefig(p, PATH_FIGS*"Fig_3.pdf")
end

# ╔═╡ 335c7914-8e5c-488e-a60c-a3b473e7de32
let
	gr()

	dict_idx = "alpha_beta=5.0"

	id = [2, 16, 30]
	if dict_idx == "alpha_beta=10.0"
		n_max_1 = [13, 100, 119] .+ 1  # for α=10
		n_max_2 = [0, 0, 79] .+ 1
		step_ticks = [3, 20, 10]
	elseif dict_idx == "alpha_beta=5.0"
		n_max_1 = [5, 25, 49] .+ 1 
		n_max_2 = [0, 0, 0] .+ 1
		step_ticks = [1, 5, 7]
	end
	
	p_i = []
	kw = (
		seriestype=:bar,
		lw=0,
		xlabel = L"| N_{k}, N_{-k}\rangle",
		ylabel = L"P(N_k, N_{-k})",
		legend=:none,
		xrotation = 60,
		# framestyle=:box,
		titlefontsize=14,
		#fontfamily="serif",
		guidefontsize=12,
		tickfontsize=12,
		colorbar_tickfontsize=12,
		dpi=300,
		rightmargin=5Plots.mm,
		leftmargin=5Plots.mm,
	)

	
	for i in eachindex(id)
		if dict_idx == "alpha_beta=5.0" && i == length(id)
			z = [data_photon_distrib_mp[dict_idx]["p_ab"][n_max_1[i] - j, n_max_2[i] + j, id[i]] for j = 0:(n_max_1[i] - n_max_2[i])]
		else
			z = [data_photon_distrib[dict_idx]["p_ab"][n_max_1[i] - j, n_max_2[i] + j, id[i]] for j = 0:(n_max_1[i] - n_max_2[i])]
		end
		labl = round(data_photon_distrib[dict_idx]["r"][id[i]]; digits = 2)
		ticks_num = [j for j=1:step_ticks[i]:(n_max_1[i] - n_max_2[i])+1]
		ticks_let = [L"|%$(n_max_1[i] - j - 1), %$(n_max_2[i] + j - 1)\rangle" for j=0:step_ticks[i]:(n_max_1[i] - n_max_2[i])]
		push!(p_i, Plots.plot(z, #aspect_ratio = :equal, 
			xticks=(ticks_num, ticks_let),
			title=L"r = %$labl";
			kw...
		))
	end

	fig = Plots.plot(p_i[3], p_i[2], p_i[1],
		layout = (1, 3),
		size = (1000, 350),
		topmargin=5Plots.mm,
		bottom_margin=20.0*Plots.mm,
	)

	# Plots.savefig(fig, PATH_FIGS*"Fig_4_alpha5.0.pdf")

	fig
end

# ╔═╡ 5a2d963b-b646-471f-b0fd-c4cb118c4b27
md"
## Erasing
"

# ╔═╡ 4a42015e-fe26-45b0-b95a-c752f0b45054
# erasing
let
	gr()

	ind = 5

	alphas = data_g2["erasing"]["alpha"]
	rs = data_g2["erasing"]["r"]
	
	omega_0 = 1.0
	omega_1 = omega_0 * rs[ind]
	init_energy = omega_0 .* (2*alphas.^2 .+ 0.5)
	init_photons = 2*alphas.^2
	
	
	x = data_g2["erasing"]["phi"]
	y_1 = (data_g2["erasing"]["n_k"][:,:,ind]) # number of photons
	y_2 = (data_g2["erasing"]["n_mk"][:,:,ind]) # number of photons
	# y_2 = omega_1 .* ((data_g2["equal"]["rll1"]["n_k"][:,:,ind] .+ 
	# 		data_g2["equal"]["rll1"]["n_mk"][:,:,ind]) .+ 0.5) ./ init_energy # energy
	# y_3 = ((data_g2["equal"]["rll1"]["g2_k"][:,:,ind] .- 1.0) .* 
	# 		data_g2["equal"]["rll1"]["n_k"][:,:,ind].^2 .+ 
	# 		data_g2["equal"]["rll1"]["n_k"][:,:,ind]) # variance in 1 mode
	# y_4 = data_entropy["equal"]["entropy_ent_vn"][1,1,:]

	l = @layout([[a b c]; d])

	p_1 = plot(x, [y_1', y_2'], 
		lw=2, m=4,
		label=["n_k, 0.5" "n_k, 1.25" "n_k, 2" "n_-k, 0.5" "n_-k, 1.25" "n_-k, 2"],
		yscale=:log10
	)

end

# ╔═╡ 6a228dc1-e68c-47fc-afa0-865081025ff3
# photons left in erasing regime, compare with vacuum photons
let
	phi_ind = 5  # 1,..,9 [0, 2π]
	x = data_g2["erasing"]["r"]
	y_1 = data_g2["erasing"]["n_mk"][:, phi_ind,:]
	y_2 = n_photons_vacuum.(u.(x), v.(x)) / 2
	
	plot(x, [y_1' y_2])
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PGFPlotsX = "8314cec4-20b6-5062-9cdb-752b83310925"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
BenchmarkTools = "~1.3.2"
FileIO = "~1.16.1"
HDF5 = "~0.17.1"
LaTeXStrings = "~1.3.1"
PGFPlotsX = "~1.6.0"
Plots = "~1.39.0"
PlutoUI = "~0.7.54"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "30fc3ce832f4486347619e142ce0e1f142c4791d"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "8cfa272e8bdedfa88b6aefbbca7c19f1befac519"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefaultApplication]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "c0dfa5a35710a193d83f03124356eef3386688fc"
uuid = "3f0dd361-4fe0-5fc6-8523-80b14ec94d85"
version = "1.1.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "26407bd1c60129062cec9da63dc7d08251544d53"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.1"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "38c8874692d48d5440d5752d6c74b0c6b0b60739"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.2+1"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8ecb0b34472a3c98f945e3c75fc7d5428d165511"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.9.3+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "8a5b4d2220377d1ece13f49438d71ad20cf1ba83"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.1.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "8f6af051b9e8ec597fa09d8885ed79fd582f33c9"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.10"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "6979eccb6a9edbbb62681e158443e79ecc0d056a"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.3.1+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b01beb91d20b0d1312a9471a36017b5b339d26de"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "PMIx_jll", "TOML", "Zlib_jll", "libevent_jll", "prrte_jll"]
git-tree-sha1 = "694458ae803b684f09c07f90459cb79655fb377d"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.0+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc6e1927ac521b659af340e0ca45828a3ffc748f"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.12+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PGFPlotsX]]
deps = ["ArgCheck", "Dates", "DefaultApplication", "DocStringExtensions", "MacroTools", "OrderedCollections", "Parameters", "Requires", "Tables"]
git-tree-sha1 = "3e7a0345b9f37da2cd770a5d47bb5cb6e62c7a81"
uuid = "8314cec4-20b6-5062-9cdb-752b83310925"
version = "1.6.0"

    [deps.PGFPlotsX.extensions]
    ColorsExt = "Colors"
    ContourExt = "Contour"
    MeasurementsExt = "Measurements"
    StatsBaseExt = "StatsBase"

    [deps.PGFPlotsX.weakdeps]
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    Contour = "d38c429a-6771-53c6-b99e-75d170b6e991"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[[deps.PMIx_jll]]
deps = ["Artifacts", "Hwloc_jll", "JLLWrappers", "Libdl", "Zlib_jll", "libevent_jll"]
git-tree-sha1 = "8b3b19351fa24791f94d7ae85faf845ca1362541"
uuid = "32165bc3-0280-59bc-8c0b-c33b6203efab"
version = "4.2.7+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a935806434c9d4c506ba941871b327b96d41f2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.0"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "1fbeaaca45801b4ba17c251dd8603ef24801dd84"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.2"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "242982d62ff0d1671e9029b52743062739255c7e"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.18.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "da69178aacc095066bad1f69d2f59a60a1dd8ad1"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.0+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522b8414d40c4cbbab8dee346ac3a09f9768f25d"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.5+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eddd19a8dea6b139ea97bdc8a0e2667d4b661720"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.0.6+1"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libevent_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll"]
git-tree-sha1 = "f04ec6d9a186115fb38f858f05c0c4e1b7fc9dcb"
uuid = "1080aeaf-3a6a-583e-a51c-c537b09f60ec"
version = "2.1.13+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.prrte_jll]]
deps = ["Artifacts", "Hwloc_jll", "JLLWrappers", "Libdl", "PMIx_jll", "libevent_jll"]
git-tree-sha1 = "5adb2d7a18a30280feb66cad6f1a1dfdca2dc7b0"
uuid = "eb928a42-fffd-568d-ab9c-3f5d54fc65b9"
version = "3.0.2+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─12f49f28-7d82-11ee-18e0-8f83344ebab3
# ╟─c5bdf327-a439-4b3e-95c0-de3f8d45c970
# ╟─9a8f509a-e8ce-4324-9c78-7fe13b4d3af1
# ╟─57b440e0-3324-4a03-b2a1-e1b8cdfb986a
# ╟─c5f0ba7a-8b39-4f26-8579-232f95a05279
# ╠═766af561-491d-40d9-92db-80fd2956bfa4
# ╠═a7eaa2b8-0b68-4987-8d3f-334595cf1891
# ╠═52cef578-5405-49a1-853d-58297fd15e25
# ╠═50ae134f-133a-4f65-8cf4-ff2256300d54
# ╠═63078392-10a0-458d-b80b-d6fc3e34dfcf
# ╠═308bc872-f072-4399-815b-3a2fb562a2d4
# ╟─3147dee6-4685-42a8-89c5-79028e8cc3a6
# ╠═0f9f6870-b5a7-4461-b58f-071889d7bc33
# ╠═2a25376b-ee9f-4172-80ae-fba2c744f70f
# ╠═15b0c4fc-7765-4036-b42a-982ec8886cc9
# ╟─79439794-d465-4ec8-ab81-27ef04110f51
# ╟─8320039b-776b-493d-a2d0-4ee56dc0d7a5
# ╠═08c2d5fe-e42c-416a-8eb5-90ed956056d3
# ╠═1c1ce53f-03c8-4ae3-b40f-6060b558e953
# ╟─5c157c15-ad97-4a1a-bec5-81975118a8ab
# ╟─9212af44-2ce8-4459-bab2-b2aabfc2b1a4
# ╠═9c70ceef-d560-4d2b-b791-27dae6eada90
# ╟─44500e2a-4eb7-401a-bae0-b2008b7751f2
# ╠═1e8a8993-3e0a-4b21-b2b2-8fa85d8994a3
# ╠═e30ab803-024c-47b4-9936-d2dfa82e4dae
# ╟─0b0f4ffc-f83f-47e0-8bb6-c99ccdc256ab
# ╟─47c3e694-9e44-4adb-9170-1b18ca10a2f4
# ╠═e065559e-d2ff-4734-a8e0-11d69c58f01a
# ╠═3da71126-bc76-43e9-9004-37a8f37d74d5
# ╠═ef6bb4ad-5ef4-47ef-ba88-f07fae9c3a1a
# ╠═335c7914-8e5c-488e-a60c-a3b473e7de32
# ╟─5a2d963b-b646-471f-b0fd-c4cb118c4b27
# ╠═4a42015e-fe26-45b0-b95a-c752f0b45054
# ╠═6a228dc1-e68c-47fc-afa0-865081025ff3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
