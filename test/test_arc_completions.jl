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

srg_completed_left = add_left_arcs(srg_upperhalf, N)

srg_completed_right = add_right_arcs(srg_upperhalf, N)

plt = plot()
plot_srg!(srg_upperhalf)
plot!(srg_completed_left)
plot!(srg_completed_right)
xlims!((-3,3))
ylims!((-3,3))
display(plt)

## now test optimal arc completions in products

srg1 = compute_lti_srg_boundary_upperhalfplane(G1, alphas, phis, frequencies)
srg2 = compute_lti_srg_boundary_upperhalfplane(G2, alphas, phis, frequencies)

N, n_bins = 20, 100
inner, outer = optimal_arc_product(srg1, srg2, N, n_bins)

srg1l = add_left_arcs(srg1, N)
srg1r = add_right_arcs(srg1, N)
srg2l = add_left_arcs(srg2, N)
srg2r = add_right_arcs(srg2, N)

plot()
plot!([inner; reverse(outer)])
plot_srg!(srg1)
plot_srg!(srg2)
plot!(srg2l)
plot!(srg2r)

# sanity check - should be the same value
maximum(abs.(srg1))*maximum(abs.(srg2))
maximum(abs.([inner; reverse(outer)]))