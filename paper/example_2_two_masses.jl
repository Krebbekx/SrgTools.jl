# Example of two masses with three nonlinear springs

using SrgTools, Plots, LaTeXStrings, ControlSystemsBase

# MSD parameters

k1trans, k2trans = 0.5, 0 # shifted centered circle
# k1trans, k2trans = 0.5-.5, 0-0.5 # left of -1 in lhp

m1, m2 = 0.5, 3
k1, k2 = 1+k1trans, 2+k1trans
d1, d2 = 0.3, 1
d12 = 1;
k12 = 0.5+2*k2trans;

s = tf("s") # define the laplace variable

A = [0 1 0 0; 
     (-k1-k12)/m1 (-d1-d12)/m1 k12/m1 d12/m1; 
     0 0 0 1; 
     k12/m2 d12/m2 (-k2-k12)/m2 (-d2-d12)/m2];

B = [0 0 0 0 0; 
     1/m1 0 1/m1 1/m1 0; 
     0 0 0 0 0; 
     0 1/m2 -1/m2 0 1/m2];

C = [1 0 0 0; 
     0 0 1 0; 
     1 0 -1 0; 
     1 0 0 0;
     0 0 1 0]

sys = ss(A,B,C,0) # define state-space system
Glfr = tf(sys)    # convert to tf

@assert(maximum(real(pole(Glfr)))<0) # check stability

Gzw = Glfr[1:3,1:3]
Gzu = Glfr[1:3,4:5]
Gyw = Glfr[4:5,1:3]
Gyu = Glfr[4:5,4:5]

# apply loop transformation scaling 
Gzw = Gzw*[1 0 0; 0 1 0; 0 0 2] 
Gyw = Gyw*[1 0 0; 0 1 0; 0 0 2]

frequencies = exp10.(collect(-5:0.01:5)) # logspace frequencies
alphas = collect(-5:.01:5)
phis = collect(0:0.005:pi)

# build circle for the SRG of the static NL
angles_circle = collect(range(0, stop=2*pi, length=300))
mu, lambda = -1/2, 1/2
radius_circle = (lambda-mu)/2 
center_circle = (lambda+mu)/2
circle = center_circle .+ radius_circle .* exp.(1im*angles_circle)
circle_inv = invert_srg(circle)

# trick to plot the complex plane except for the inverted circle as the SRG of the inverted circle
farfield = 5
circle_inv_plot = [circle_inv; farfield; farfield+farfield*1im; -farfield+farfield*1im; -farfield-farfield*1im; farfield-farfield*1im; farfield]

plt = plot()
plot_srg!(circle_inv_plot)
plot_lti_srg!(Gzw, alphas, phis, frequencies)
xlims!((-3,3))
ylims!((-2.5,2.5))
annotate!(-1.7, 2, text(L"$\operatorname{SRG}(\tilde{\Phi})^{-1}$", 22))
annotate!(.5, -0.9, text(L"$\operatorname{SRG}(\tilde{G}_\mathrm{zw})$", 22))
savefig("paper/figures/two_masses_Gzw_Phi.pdf") 
display(plt)

# compute SRGs for each individual operator
plot()
plot_lti_srg!(Gyu, alphas, phis, frequencies)
savefig("paper/figures/example_2_srg_Gyu.pdf") 

plot()
plot_lti_srg!(Gzu, alphas, phis, frequencies)
savefig("paper/figures/example_2_srg_Gzu.pdf") 

plot()
plot_lti_srg!(Gyw, alphas, phis, frequencies)
savefig("paper/figures/example_2_srg_Gyw.pdf") 


# compute the srg of (phi^{-1} - Gzw) in case of inner and outer radius
srg_Gzw_upperhalf = compute_lti_srg_boundary_upperhalfplane(Gzw, alphas, phis, frequencies)
circle_inv_upperhalf = [z for z in circle_inv if imag(z) >= 0]

# note that full srg, not only upperhalfplane, is needed for sum operation
srg_Gzw_chord = add_chords(srg_Gzw_upperhalf, 10) # add chords
srg_denom = [z1+z2 for z1 in circle_inv for z2 in -srg_Gzw_chord] # minus sign only with real numbers

# compute inner/outer radius of upper half plane part of denom. Only use inner radius
inner_rad_denom, _ = inner_outer_radius_srg_binning([z for z in srg_denom if imag(z) >= 0], 100);

srg_return_upperhalf = invert_srg(inner_rad_denom)
srg_return = [srg_return_upperhalf; reverse(conj(srg_return_upperhalf))]


plot()
plot_srg!(srg_return)

srg_Gzu = compute_lti_srg_boundary_upperhalfplane(Gzu, alphas, phis, frequencies)
srg_Gyw = compute_lti_srg_boundary_upperhalfplane(Gyw, alphas, phis, frequencies)
srg_Gyu = compute_lti_srg_boundary_upperhalfplane(Gyu, alphas, phis, frequencies)
N, n_bins = 20, 100

inner1, outer1 = improved_arc_product(srg_return_upperhalf, srg_Gzu, N, n_bins);
srg1 = [inner1; outer1];
inner2, outer2 = improved_arc_product(srg_Gyw, srg1, N, n_bins);
srg2 = [inner2; outer2];
inner, outer = improved_chord_sum(srg2, srg_Gyu, N, n_bins);

srg_upperhalf = [inner; reverse(outer)];
maximum(abs.(srg_upperhalf)) # gives the closed-loop incremental L2-gain bound

plot()
plot_srg!([outer; reverse(conj(outer))])
savefig("paper/figures/example_2_srg_full_lfr.pdf") 

# compute rough gain bound using max gain of each operator, quicker but more conservative
maximum(abs.(srg_Gyu)) + maximum(abs.(srg_Gzu))*maximum(abs.(srg_Gyw))*maximum(abs.(srg_return))  # worst case bound


### Alternatively, build tf directly
# abbreviate G_{y_i u_j} = Gyiuj
# denominator = (k12+d12*s)^2-(k1+k12+s*(d1+d12+m1*s))*(k12+k2+s*(d12+d2+m2*s))

# @assert(maximum(real(pole(1/denominator)))<0) # check stability

# Gy1u1 = -(k12 + k2 + s*(d12 + d2 + m2*s))/denominator
# Gy1u2 = -(k12 + d12*s)/denominator
# Gy2u1 = -(k12 + d12*s)/denominator
# Gy2u2 = -(k1 + k12 + s*(d1 + d12 + m1*s))/denominator

# Gz1w3 = Gy1w3 = Gy1u1-Gy1u2 # pos x1 from NL k12 spring force
# Gz2w3 = Gy2w3 = Gy2u1-Gy2u2 # pos x2 from NL k12 spring force
# Gz3u1 = Gz3w1 = Gy1u1-Gy2u1 # pos diff. x1-x2 from force on m1
# Gz3u2 = Gz3w2 = Gy1u2-Gy2u2 # pos diff. x1-x2 from force on m2

# Gz3w3 = Gz1w3-Gz2w3         # pos diff. x1-x2 from NL k12 spring force

# Glfr = [Gy1u1 Gy1u2 Gz1w3 Gy1u1 Gy1u2;
#         Gy2u1 Gy2u2 Gz2w3 Gy2u1 Gy2u2;
#         Gz3w1 Gz3w2 Gz3w3 Gz3u1 Gz3u2;
#         Gy1u1 Gy1u2 Gy1w3 Gy1u1 Gy1u2;
#         Gy2u1 Gy2u2 Gy2w3 Gy2u1 Gy2u2]