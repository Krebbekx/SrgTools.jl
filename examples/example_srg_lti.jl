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

# this script shows how to plot the SRG of an LTI system G, and to obtain the points that describe the boundary of the SRG as a polygon

using SrgTools, Plots


G(s::Number) = [s/(s+1)     2*s/(s+1); 
                s^2/(s+1)^2 s/(s+1)^2] # 2x2 transfer function

frequencies = exp10.(collect(-3:0.01:3)) # grid of frequency points to evaluate over

alphas = collect(-5:.05:5) # points on the real axis points with respect to which the SRG radius is computed

phis = collect(0:0.01:pi) # angular resolution of the outer bound of the SRG in the upper half plane. In this case, it is the upper half plane since [0,pi] is the interval, with resolution of 0.01 rad.

# animation of the max (blue) and min (blue) gain circles
plt = plot()
plot_lti_srg_circles!(G, alphas, phis, frequencies)

# plot of the resulting SRG
plt = plot()
plot_lti_srg!(G, alphas[1:end], phis, frequencies)
xlims!((-3,3))
ylims!((-3,3))
display(plt)

savefig("examples/myplot.pdf")   # save the result

# to obtain the complex points that parameterize the boundary (as a polygon) of the SRG in the upper half complex plane
srg_upperhalf = compute_lti_srg_boundary_upperhalfplane(G, alphas, phis, frequencies)