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

# Helper functions for frequently occurring computations in the complex plane.


"""
Given two polar decompositions of the same number z, defined by 

    z = z0 + r0 exp{j phi} = alpha + r_alpha exp{j varphi}, (1)

such that |z0-alpha| < r_alpha and phi in [0, pi].

This function computes varphi using a direct computation with quadratic polynomials. The reason to use this function is computational speed: it eliminates the need to solve Eq. (1) above. 
"""
function compute_varphi_srg_center(z0, alpha, r_alpha, phi)
    # test if computation is allowed
    @assert(0 <= phi <= pi)

    # special case
    if phi == pi/2
        return acos(-(alpha-z0)/r_alpha)
    elseif phi > pi/2 # go to the right atan2 sector
        phi = phi-pi
    end

    # fill in abc formula
    tanphi2 = tan(phi)^2
    a = r_alpha^2*tanphi2 + r_alpha^2;
    b = 2*r_alpha*(alpha - z0)*tanphi2;
    c = (alpha - z0)^2*tanphi2 - r_alpha^2;

    x1, x2 = (-b-sqrt(b^2-4*a*c))/(2*a), (-b+sqrt(b^2-4*a*c))/(2*a)

    @assert(abs(x1)<=1 || abs(x2)<=1)

    if abs(x1) > 1
        return acos(x2)
    elseif  abs(x2) > 1
        return acos(x1)
    else
        # two possible solutions due to equation squaring. Time to pick the right one
        errorfun(x) = tan(phi)*(alpha - z0 + r_alpha*x) - r_alpha*sqrt(1 - x^2)
        if abs(errorfun(x1)) < abs(errorfun(x2))
            return acos(x1)
        else
            return acos(x2)
        end
    end
end


"""
Computes the inner circle (closest points to origin) of an SRG traced along a circle with radius r and center c with N points.
"""
function srg_plus_circle_boundary(srg::Vector, r, c, N)
    phis = collect(range(0, stop=2*pi, length=N))

    circle = c .+ r .* exp.(1im*phis)

    boundary = zeros(Complex, N)

    for i in eachindex(circle)
        mindex = argmin(abs.(circle[i] .+ srg)) # compute point closest to the origin
        boundary[i] = circle[i] .+ srg[mindex]
    end 

    return boundary
end


"""
Sorts the complex points in n_bins evenly spaced bins of arguments between the angles start and stop in radians.
"""
function sort_points_in_bins(points::Vector, n_bins::Int, start::Number, stop::Number)
    indices_sorted = [[] for j=1:n_bins]        # keep track of which point is in which bin
    
    phis = collect(range(start, stop=stop, length = n_bins+1))
    args = imag(log.(points)) # precompute all arguments

    for i in eachindex(points)                         # first, sort all numbers in the bins
        for j = 1:n_bins
            if phis[j] <= args[i] < phis[j+1]
                push!(indices_sorted[j], i)     # add index to bin, break out of loop
                break
            end 
        end 
    end

    return indices_sorted
end


"""
Computes the inner and outer radius (w.r.t. origin) as sets of complex points that bound the SRG away from the origin.
"""
function inner_outer_radius_srg_binning(srg::Vector, n_bins::Int)
    @assert(n_bins >= 1)

    srg_arg = imag(log.(srg))                       # compute arguments

    start, stop = minimum(srg_arg), maximum(srg_arg)

    srg_indices_sorted = sort_points_in_bins(srg, n_bins, start, stop) # sort srg points in bins or arguments

    # compute the indices and magnitudes in the non-empty bins
    inner_radius_indices = []
    outer_radius_indices = []

    # compute the min and max index in each non-empty
    for j = 1:n_bins
        if length(srg_indices_sorted[j]) > 0
            absvals = abs.(srg[srg_indices_sorted[j]]) # compute magnitudes in the bin
            push!(inner_radius_indices, srg_indices_sorted[j][argmin(absvals)]) 
            push!(outer_radius_indices, srg_indices_sorted[j][argmax(absvals)]) 
        end 
    end 

    # fetch the points corresponding to these indices
    inner_radius, outer_radius = srg[inner_radius_indices], srg[outer_radius_indices] 

    return inner_radius, outer_radius
end