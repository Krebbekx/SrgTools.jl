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

# Example of two two MIMO transfer functions

using SrgTools, Plots, LaTeXStrings, ControlSystemsBase


s = tf("s") # define the laplace variable

G1 = [s/(s+1) s^2/(s^2+s+1) 1/(2*s+1)]

G2 = [s^2/(s^2+s+1) 1/(2*s+1); (s+1)/((s+3)*(s^2+s+1)) (s+3)/(s+1); (s^2-1)/((s+3)*(s+2)) s/(s+2)]

frequencies = exp10.(collect(-5:0.01:5)) # logspace frequencies
alphas = collect(-5:.01:5)
phis = collect(0:0.005:pi)

plt = plot()
plot_lti_srg!(G1, alphas, phis, frequencies)
# xlims!((-1.5,2))
# ylims!((-1.5,1.5))
savefig("paper/figures/example_lti_G1.pdf") 
display(plt)

plt = plot()
plot_lti_srg!(G2, alphas, phis, frequencies)
# xlims!((-3,3))
# ylims!((-2.5,2.5))
savefig("paper/figures/example_lti_G2.pdf") 
display(plt)