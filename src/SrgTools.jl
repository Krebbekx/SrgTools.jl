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

module SrgTools

#=
    This module contains functions to manipulate SRGs as vectors of complex numbers. functionalities include
    - computing SRGs for LTI systems
    - computing the SRG of inverses, sums and products of SRGs
    - plotting SRGs    

    The intended purpose of this toolbox is to use it as a "sandbox" for SRG computations. The user has to verify the SRG separation by hand, and compute separations or SRG radii using the native Julia functionality: maximum(abs.(srg)) 
    
    Topic of future work is to develop these tools further such that they can serve as the cost function in optimization algorithms (e.g. controller/system as input, gain bound obtained via SRG computations as output). 
=#


using Plots, LinearAlgebra, ControlSystemsBase


export 
    plot_srg!,
    invert_srg,
    compute_lti_srg_boundary_upperhalfplane,
    plot_lti_srg!,
    inner_outer_radius_srg_binning,
    add_chords,
    add_left_arcs,
    add_right_arcs,
    improved_chord_sum,
    improved_arc_product,
    plot_lti_srg_circles!


include("complex_geometry.jl")      # functions for common computations with complex numbers
include("lti_tools.jl")             # functions to compute the SRG of LTI systems
include("interconnections.jl")      # functions to numerically handle sums, products and inverses of SRGs
include("plotting.jl")              # all plotting functionalities


end # module SrgTools
