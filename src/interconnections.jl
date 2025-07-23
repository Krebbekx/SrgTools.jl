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

# this file contains functions to numerically evaluate SRG interconnections: inverses, sums and products. Also, functions to add chords and arcs, possibly in an improved fashion.

"""
Function to invert the SRG using a Möbius inverse. 

# Arguments
- `srg::Vector`: points in the complex plane that parameterize the boundary of the SRG as polygon vertices.
"""
invert_srg(srg::Vector) = conj(1 ./ srg)


"""
Adds chords to the srg. Chords consist of N points. The result is an unsorted point cloud.
"""
function add_chords(srg::Vector, N::Int)
    srg_upperhalf = [z for z in srg if imag(z) >= 0] # filter out lower half plane

    srg_completed = ComplexF64[]

    for z in srg_upperhalf
        if imag(z) > 0
            alphas = collect(range(0, stop=1, length=N))
            chord = [z + (conj(z)-z)*alpha for alpha in alphas]
            append!(srg_completed, chord)
        else
            push!(srg_completed, z) # just add the number
        end 
    end 

    return srg_completed
end


"""
Adds left-hand arcs to the srg. Arcs consist of N points. The result is an unsorted point cloud.
"""
function add_left_arcs(srg::Vector, N::Int)
    srg_upperhalf = [z for z in srg if imag(z) >= 0] # filter out lower half plane

    srg_completed = ComplexF64[]

    for z in srg_upperhalf
        if imag(z) > 0 || real(z) > 0
            phi = imag(log(z))      # compute angle
            phis = collect(range(phi, stop=2*pi-phi, length=N))
            r = abs(z)              # compute magnitude
            arc = r*exp.(1im*phis)
            append!(srg_completed, arc)
        else
            push!(srg_completed, z) # just add the number
        end 
    end 

    return srg_completed
end


"""
Adds right-hand arcs to the srg. Arcs consist of N points. The result is an unsorted point cloud.
"""
function add_right_arcs(srg::Vector, N::Int)
    srg_upperhalf = [z for z in srg if imag(z) >= 0] # filter out lower half plane

    srg_completed = ComplexF64[]

    for z in srg_upperhalf
        if imag(z) > 0 || real(z) < 0
            phi = imag(log(z))      # compute angle
            phis = collect(range(-phi, stop=phi, length=N))
            r = abs(z)              # compute magnitude
            arc = r*exp.(1im*phis)
            append!(srg_completed, arc)
        else
            push!(srg_completed, z) # just add the number
        end 
    end 

    return srg_completed
end


"""
Computes SRG sum by computing the two chord completions and outputting the largest minradius and smallest maxradius.
- N      : amount of points per chord
- n_bins : amount of bins

srg1 and srg2 may be upperhalfplane only.

Outputs improved chord completed upperhalfplane.
"""
function improved_chord_sum(srg1::Vector, srg2::Vector, N::Int, n_bins::Int)
    # compute all completions
    srg1c = add_chords(srg1, N)
    srg2c = add_chords(srg2, N)

    # compute all products
    p_1c2 = [z1+z2 for z1 in srg1c for z2 in srg2]
    p_12c = [z1+z2 for z1 in srg1 for z2 in srg2c]

    # flip back to upper half plane
    p_1c2 = [z for z in [p_1c2; conj(p_1c2)] if imag(z) >= 0]
    p_12c = [z for z in [p_12c; conj(p_12c)] if imag(z) >= 0]

    points = [p_1c2, p_12c] # combine all points in one array

    # sort points in upper half plane
    i_1c2 = sort_points_in_bins(p_1c2, n_bins, 0, pi)
    i_12c = sort_points_in_bins(p_12c, n_bins, 0, pi)

    inner_radius = ComplexF64[]
    outer_radius = ComplexF64[]

    # loop over bins
    for i = 1:n_bins
        # for each completion, bins must be non-empty, otherwise the non-empty one wins
        if length(i_1c2[i])>0  && length(i_12c[i])>0 
            # compute for each of the two sets the index of the max and min point
            inner_indices = [i_1c2[i][argmin(abs.(p_1c2[i_1c2[i]]))], i_12c[i][argmin(abs.(p_12c[i_12c[i]]))]]
            outer_indices = [i_1c2[i][argmax(abs.(p_1c2[i_1c2[i]]))], i_12c[i][argmax(abs.(p_12c[i_12c[i]]))]]

            # compute which set of the four contains the smallest outer radius and largest inner radius
            min_index = argmax([abs(points[i][inner_indices[i]]) for i=1:2])
            max_index = argmin([abs(points[i][outer_indices[i]]) for i=1:2])

            push!(inner_radius, points[min_index][inner_indices[min_index]])
            push!(outer_radius, points[max_index][outer_indices[max_index]])
        end 
    end 

    return inner_radius, outer_radius
end


"""
Computes SRG multiplication by computing all 4 arc combinations and outputting the largest minradius and smallest maxradius.
- N      : amount of points per arc
- n_bins : amount of bins

srg1 and srg2 may be upperhalfplane only.

Outputs improved arc completed upperhalfplane.
"""
function improved_arc_product(srg1::Vector, srg2::Vector, N::Int, n_bins::Int)
    # compute all completions
    srg1_left = add_left_arcs(srg1, N)
    srg1_right = add_right_arcs(srg1, N)
    srg2_left = add_left_arcs(srg2, N)
    srg2_right = add_right_arcs(srg2, N)

    # compute all products
    p_1l2 = [z1*z2 for z1 in srg1_left for z2 in srg2]
    p_1r2 = [z1*z2 for z1 in srg1_right for z2 in srg2]
    p_12l = [z1*z2 for z1 in srg1 for z2 in srg2_left]
    p_12r = [z1*z2 for z1 in srg1 for z2 in srg2_right]

    # flip back to upper half plane
    p_1l2 = [z for z in [p_1l2; conj(p_1l2)] if imag(z) >= 0]
    p_1r2 = [z for z in [p_1r2; conj(p_1r2)] if imag(z) >= 0]
    p_12l = [z for z in [p_12l; conj(p_12l)] if imag(z) >= 0]
    p_12r = [z for z in [p_12r; conj(p_12r)] if imag(z) >= 0]

    points = [p_1l2, p_1r2, p_12l, p_12r] # combine all points in one array

    # sort points in upper half plane
    i_1l2 = sort_points_in_bins(p_1l2, n_bins, 0, pi)
    i_1r2 = sort_points_in_bins(p_1r2, n_bins, 0, pi)
    i_12l = sort_points_in_bins(p_12l, n_bins, 0, pi)
    i_12r = sort_points_in_bins(p_12r, n_bins, 0, pi)

    inner_radius = ComplexF64[]
    outer_radius = ComplexF64[]

    # loop over bins
    for i = 1:n_bins
        # for each completion, bins must be non-empty, otherwise the non-empty one wins
        if length(i_1l2[i])>0 && length(i_1r2[i])>0 && length(i_12l[i])>0 && length(i_12r[i])>0
            # compute for each of the four sets the index of the max and min point
            inner_indices = [i_1l2[i][argmin(abs.(p_1l2[i_1l2[i]]))], i_1r2[i][argmin(abs.(p_1r2[i_1r2[i]]))], i_12l[i][argmin(abs.(p_12l[i_12l[i]]))], i_12r[i][argmin(abs.(p_12r[i_12r[i]]))]]
            outer_indices = [i_1l2[i][argmax(abs.(p_1l2[i_1l2[i]]))], i_1r2[i][argmax(abs.(p_1r2[i_1r2[i]]))], i_12l[i][argmax(abs.(p_12l[i_12l[i]]))], i_12r[i][argmax(abs.(p_12r[i_12r[i]]))]]

            # compute which set of the four contains the smallest outer radius and largest inner radius
            min_index = argmax([abs(points[i][inner_indices[i]]) for i=1:4])
            max_index = argmin([abs(points[i][outer_indices[i]]) for i=1:4])

            push!(inner_radius, points[min_index][inner_indices[min_index]])
            push!(outer_radius, points[max_index][outer_indices[max_index]])
        end 
    end 

    return inner_radius, outer_radius
end