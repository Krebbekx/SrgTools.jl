# =====================================================================
# This file is part of the SrgTools.jl project, a software toolbox 
# written in Julia for computations with (MIMO) Scaled Relative Graphs.
#
# Copyright (c) 2025, Julius P.J. Krebbekx, Eindhoven, The Netherlands
# All rights reserved.
#
# This source code is licensed under the BSD 3-Clause License.
# See the LICENSE file in the project root for license information.
# =====================================================================

using SrgTools, Plots

G1(s::Number) = [s/(s+1) 2*s/(s+1); s^2/(s+1)^2 s/(s+1)^2] # 2x2 transfer function

G2(s::Number) = [s/(s+1) -s/(s+1); s/(s+1)^2 s/(s+1)] # 2x2 transfer function

frequencies = exp10.(collect(-3:0.01:3)) # grid of frequency points to evaluate over

alphas = collect(-5:.05:5) # real axis points to compute the SRG radius with 

phis = collect(0:0.01:pi) # angular resolution of the outer bound of the SRG in the upper half plane

# to obtain the complex points that parameterize the boundary (as a polygon) of the SRG in the upper half complex plane
srg_upperhalf = compute_lti_srg_boundary_upperhalfplane(G1, alphas, phis, frequencies)

N = 100

srg_completed = add_chords(srg_upperhalf, N)

plt = plot()
plot_srg!(srg_upperhalf)
plot!(srg_completed)
xlims!((-3,3))
ylims!((-3,3))
display(plt)

# now test optimal chord completions

srg1 = compute_lti_srg_boundary_upperhalfplane(G1, alphas, phis, frequencies)
srg2 = compute_lti_srg_boundary_upperhalfplane(G2, alphas, phis, frequencies)

# other example to show difference!
# srg2 = [1+0.1*exp(1im*phi) for phi in range(0, stop=2*pi, length=200)]

N, n_bins = 10, 100
inner_opt, outer_opt = optimal_chord_sum(srg1, srg2, N, n_bins)

srg1c = add_chords(srg1, N)
srg2c = add_chords(srg2, N)

inner, outer = inner_outer_radius_srg_binning([z1+z2 for z1 in srg1c for z2 in srg2c], n_bins);

plot()
plot!([inner_opt; reverse(outer_opt)])
plot_srg!(srg1)
plot_srg!(srg2)
plot!(srg1c)
plot!(srg2c)
plot!([inner; reverse(outer)])

# sanity check - should be the same value
maximum(abs.(srg1))*maximum(abs.(srg2))
maximum(abs.([inner; reverse(outer)]))