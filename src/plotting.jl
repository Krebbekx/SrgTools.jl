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

# Contains all SRG plotting functionalities.

"""
Function to plot SRGs as a set in the complex plane. Must act on an existing plot.

# Arguments
- `srg::Vector`: points in the complex plane that parameterize the boundary of the SRG as polygon vertices.

Optional arguments for fillcolor and fillalpha with default values :black and 0.4, respectively.
"""
function plot_srg!(srg::Vector, fillcolor=:black, fillalpha=0.4)
    x, y = real(srg), imag(srg)
    # plot!(Shape(x, y), fillcolor = RGB(0/255, 75/255, 200/255), fillalpha = 0.7, linealpha = 0, grid = false, aspect_ratio = :equal, legend = false) # other color choice

    # plot!(Shape(x, y), fillcolor = :black, fillalpha = 0.4, linealpha = 0, grid = false, aspect_ratio = :equal, legend = false)

    plot!(Shape(x, y), fillcolor = fillcolor, fillalpha = fillalpha, linealpha = 0, grid = false, aspect_ratio = :equal, legend = false)


    xlabel!("Re")
    ylabel!("Im")
    plot!(  xlabel = "Re",
            label = "Im",
            xguidefont = font(20, "Computer Modern"),   # Axis label font size and family
            yguidefont = font(20, "Computer Modern"),
            xtickfont = font(14, "Computer Modern"),    # Tick label font size and family
            ytickfont = font(14, "Computer Modern"),
            grid = true
            )
end 


"""
Plots intersection of max gain circles and plots white min gain circles on top.

Removes the circles with centers outside of [a, b]

G can be either a scalar/vector/matrix function or a TransferFunction object from ControlSystemsBase.

Optional arguments for fillcolor and fillalpha with default values :black and 0.4, respectively.
"""
function plot_lti_srg!(G::Union{Function, TransferFunction}, alphas::Vector, phis::Vector, frequencies::Vector, fillcolor=:black, fillalpha=0.4)
    srg_upperhalf = compute_lti_srg_boundary_upperhalfplane(G, alphas, phis, frequencies)
    srg = [srg_upperhalf; conj(reverse(srg_upperhalf))] # build full SRG

    plot_srg!(srg, fillcolor, fillalpha)
end


function plot_lti_srg_circles!(phis::Vector, alphas::Vector, max_radii::Vector)
    n = length(alphas)
    for i = 1:n
        z = [alphas[i] .+ max_radii[i]*exp.(1im*phi) for phi in phis]
        display(plot!(real(z), imag(z)))
    end 
end 


"""
OLD 

Plots intersection of max gain circles and plots white min gain circles on top.
"""
function plot_lti_srg_circles!(G::Union{Function, TransferFunction}, alphas::Vector, phis::Vector, frequencies::Vector)
    min_radii, max_radii = compute_srg_circles_lti(G, alphas, frequencies)

    # generate a point z0 that is contained in all max_radii circles
    a, b = find_interval(alphas, max_radii)
    z0 = (b+a)/2
    print(a)
    print("\n")
    print(b)
    print("\n")
    print(z0)

    srg_max_radius = boundary_max_gain_circle_intersections(z0, alphas, max_radii, phis)

    srg_upperhalf = z0.+srg_max_radius.*exp.(1im*phis)  
    srg = [srg_upperhalf; conj(reverse(srg_upperhalf))]
    # display(plot!(Shape(real(srg), imag(srg))))
    display(plot!(Shape(real(srg), imag(srg)), fillcolor = :black, fillalpha = 0.2, linealpha = 0.0, grid = false, aspect_ratio = :equal, legend = false))

    xlabel!("Re")
    ylabel!("Im")

    circle(base, radius) = [base .+ radius*exp.(1im*angle) for angle in collect(range(0, stop = 2*pi, length = 100))]

    for i in eachindex(alphas)
        z = circle(alphas[i], min_radii[i]) # create circle with base alpha and radius r_alpha
        display(plot!(Shape(real(z), imag(z)), fillalpha=0, linecolor = :black)) # plot a white circle
        z = circle(alphas[i], max_radii[i]) # create circle with base alpha and radius r_alpha
        display(plot!(Shape(real(z), imag(z)), fillalpha=0, linecolor = :blue)) # plot a white circle
    end 
end