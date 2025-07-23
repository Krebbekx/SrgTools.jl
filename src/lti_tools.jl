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

# Contains all functions to compute the SRG of an LTI function. 

"""
Computes the minimal and maximal singular values of the transfer matrix G over a set of frequency points in frequencies.

G can be either a scalar/vector/matrix function, or a TransferFunction object from ControlSystemsBase.
"""
function compute_min_max_gain(G::Union{Function, TransferFunction}, frequencies::Vector)
    # start with the first frequency value
    sigmas = svd(G(frequencies[1]*1im)).S
    σ_min, σ_max = minimum(sigmas), maximum(sigmas)

    for omega in frequencies[2:end]
        sigmas = svd(G(omega*1im)).S

        # update largest gain if necessary
        if σ_max < maximum(sigmas)
            σ_max = maximum(sigmas)
        end

        # update smallest gain if necessary
        if σ_min > minimum(sigmas)
            σ_min = minimum(sigmas)
        end
    end

    return σ_min, σ_max
end 


"""
Function to compute the max and min radius of the SRG w.r.t. real-valued base points α on the real line.

G can be either a scalar/vector/matrix function or a TransferFunction object from ControlSystemsBase.
"""
function compute_srg_circles_lti(G::Union{Function, TransferFunction}, alphas::Vector, frequencies::Vector)

    # since G can be a function mapping to a matrix (includes row vector), vector (a qx1 dimensional matrix) or a scalar, one has to carefully cast all options to the same format before proceeding to the SVD procedure 
    if length(size(G(0))) == 2 # matrix and row vector case
        q, p = size(G(0))
    elseif length(size(G(0))) == 1 # column vector case, or scalar in []
        p = 1
        q = size(G(0))[1]
    else # scalar case
        p = q = 1
    end

    # correct for dimensions and create correct identity
    id = Matrix(I,max(q,p),p)
    G_alpha(s, alpha) = [G(s); zeros(max(p-q,0),p)] - alpha*id

    N = length(alphas) # amount of circles
    min_radii, max_radii = zeros(N), zeros(N)

    # loop over all real alpha values
    for i in 1:N
        G_a(s) = G_alpha(s, alphas[i]) # define the alpha function
        min_radii[i], max_radii[i] = compute_min_max_gain(G_a, frequencies) # compute gains
    end

    return min_radii, max_radii
end


"""
Computes for phi in phis the outer radius of the SRG, w.r.t. the real point z0 in the SRG.

It is advised to take phis as a grid in [0,pi].
"""
function boundary_max_gain_circle_intersections(z0, alphas::Vector, max_radii::Vector, phis::Vector)
    n = length(phis)
    m = length(alphas)
    srg_max_radius = zeros(n)

    radii = zeros(m) # temporary container for radii of all alpha circles
    for i = 1:n
        for j = 1:m
            varphi = compute_varphi_srg_center(z0, alphas[j], max_radii[j], phis[i])
            radii[j] = abs(alphas[j]+max_radii[j]*exp(1im*varphi)-z0)
        end
        srg_max_radius[i] = minimum(radii)
    end 

    return srg_max_radius
end


"""
Finds the interval [a,b] that is contained in all circles with center in alphas and radius in max_radii
"""
function find_interval(alphas::Vector, max_radii::Vector)
    # generate initial interval using the first circle
    a, b = alphas[1]-max_radii[1], alphas[1]+max_radii[1]

    for i = 2:length(alphas)
        a, b = max(a, alphas[i]-max_radii[i]), min(b, alphas[i]+max_radii[i])
    end 

    return a, b
end


"""
Removes all min gain circles from the max gain outer approx of the SRG.

First loops over all circles with center outside [a,b], then all remaining min gain circles.

Outputs complex contour of the SRG, which is the connected inner and outer radius w.r.t. the point z0.
"""
function remove_min_circles_from_max_gain(z0, srg_max_radius::Vector, phis::Vector, alphas::Vector, min_radii::Vector)
    srg_outer_radius = z0.+srg_max_radius.*exp.(1im*phis) # compute the max gain contour as complex numbers

    a, b = minimum(real(srg_outer_radius)), maximum(real(srg_outer_radius)) # compute min and max real values of the convex set obtained by the max_radii circles

    srg_outer_radius = [srg_outer_radius; collect(range(a, stop = b, length = length(phis)))] # return with a line on the real axis for the minimal circles

    for i in eachindex(alphas) # loop over all min_radii circles outside [a,b]
        if !(a <= alphas[i] <= b) 
            for j in eachindex(srg_outer_radius)    # per circle, loop over all points
                z = srg_outer_radius[j]-alphas[i]   # compute srg boundary point from the center of the min_radii circle
                varphi = imag(log(z))     # compute argument
                varphi = atan(imag(z), real(z))
                radius = abs(z)           # compute radius

                if radius < min_radii[i] # update value if min_radii circle intersects the SRG
                    srg_outer_radius[j] = alphas[i] + min_radii[i]*exp(1im*varphi)
                end 
            end 
        end  
    end 

    for i in eachindex(alphas) # loop over all min_radii circles inside [a,b]
        if a <= alphas[i] <= b 
            for j in eachindex(srg_outer_radius)
                z = srg_outer_radius[j]-alphas[i]
                if abs(z) < min_radii[i] # change must occur
                    re_z = real(z)
                    srg_outer_radius[j] = alphas[i] + re_z + 1im*sqrt(1-(re_z/min_radii[i])^2)*min_radii[i] # project on top of the circle
                end 
            end 
        end  
    end

    return srg_outer_radius # this is the boundary of the SRG in the upper half plane
end


"""
Computes the set of complex numbers that parameterize the boundary of the polygon that describes the upper half plane part of the SRG of an LTI operator.

G can be either a scalar/vector/matrix function or a TransferFunction object from ControlSystemsBase.
"""
function compute_lti_srg_boundary_upperhalfplane(G::Union{Function, TransferFunction}, alphas::Vector, phis::Vector, frequencies::Vector)
    min_radii, max_radii = compute_srg_circles_lti(G, alphas, frequencies)

    # generate a point z0 that is contained in all max_radii circles
    a, b = find_interval(alphas, max_radii)
    z0 = (b+a)/2

    srg_max_radius = boundary_max_gain_circle_intersections(z0, alphas, max_radii, phis)

    srg_upperhalf = remove_min_circles_from_max_gain(z0, srg_max_radius, phis, alphas, min_radii)

    return srg_upperhalf
end