using SrgTools, Plots, LaTeXStrings, ControlSystemsBase

# Nonlinear network example: controlled Lur'e system with load side nonlinearity in feedback with a linear controller.

s = tf("s") # define the laplace variable

# define the stable controller and unstable LTI plant
K = 1/(s+1)
G = 3/((s-2)*(s/10+1))

# by Nyquist criterion, k1 > 2/3 is required
k1, k2 = 2, 3 # loop transformation gains. Make sure that these result in stable plants!

Gt = G/(1+k2*G) # loop transformed G

S = 1/(1+k1*K*Gt) # sensitivity function

# elements of the LFR matrix
G_z1u = S*K
G_yu  = S*k1*K*Gt
G_z1w = -S*K*Gt
G_yw  = S*Gt

Glfr = [G_z1w G_z1w G_z1u; 
        G_yw  G_yw  G_yu ;
        G_yw  G_yw  G_yu]

Gzw = Glfr[1:2,1:2]
Gzu = Glfr[1:2,3:3]
Gyw = Glfr[3:3,1:2]
Gyu = Glfr[3:3,3:3]

frequencies = exp10.(collect(-3:0.01:3)) # logspace frequencies to evaluate the max and min singular values at
alphas = collect(-5:.02:5)               # points on the real line with respect to which the max and min gains are computed
phis = collect(0:0.01:pi)                # angular resolution of the resulting SRG

# build circle for the SRG of the static NL
angles_circle = collect(range(0, stop=2*pi, length=300))
mu, lambda = min(-k1,1-k2), max(1-k1, 2-k2)
radius_circle = (lambda-mu)/2 
center_circle = (lambda+mu)/2
circle = center_circle .+ radius_circle .* exp.(1im*angles_circle)
circle_inv = invert_srg(circle)

plt = plot()
plot_lti_srg!(Gzw, alphas, phis, frequencies)
plot_srg!(circle_inv)
xlims!((-1.5,1.3))
ylims!((-1.0,1.0))
annotate!(-.9, -0.4, text(L"$\operatorname{SRG}(\Phi)^{-1}$", 26))
annotate!(.9, -0.85, text(L"$\operatorname{SRG}(G_\mathrm{zw})$", 26))
savefig("paper/figures/example_1_k23.pdf") 
display(plt)

# compute the srg of (phi^{-1} - Gzw) in case of inner and outer radius
srg_Gzw_upperhalf = compute_lti_srg_boundary_upperhalfplane(Gzw, alphas, phis, frequencies)
circle_inv_upperhalf = [z for z in circle_inv if imag(z) >= 0]

# note that full srg, not only upperhalfplane, is needed for sum operation
srg_denom = [z1+z2 for z1 in circle_inv for z2 in -[srg_Gzw_upperhalf; conj(srg_Gzw_upperhalf)]] # minus sign only with real numbers

# compute inner/outer radius of upper half plane part of denom
inner_rad_denom, outer_rad_denom = inner_outer_radius_srg_binning([z for z in srg_denom if imag(z) >= 0], 100);

srg_return_upperhalf = invert_srg([reverse(inner_rad_denom); outer_rad_denom])
srg_return = [srg_return_upperhalf; reverse(conj(srg_return_upperhalf))]


# compute SRG of the whole LFR, including chords and arcs
srg_Gzu = compute_lti_srg_boundary_upperhalfplane(Gzu, alphas, phis, frequencies)
srg_Gyw = compute_lti_srg_boundary_upperhalfplane(Gyw, alphas, phis, frequencies)
srg_Gyu = compute_lti_srg_boundary_upperhalfplane(Gyu, alphas, phis, frequencies)
N, n_bins = 20, 100

inner1, outer1 = improved_arc_product(srg_return_upperhalf, srg_Gzu, N, n_bins);
srg1 = [inner1; outer1];
inner2, outer2 = improved_arc_product(srg_Gyw, srg1, N, n_bins);
srg2 = [inner2; outer2];
inner, outer = improved_chord_sum(srg2, srg_Gyu, N, n_bins);

srg_upperhalf = [inner; reverse(outer)]; # append the inner and outer bound of the SRG
maximum(abs.(srg_upperhalf)) # gives the closed-loop incremental L2-gain bound

plot()
plot_srg!([outer; reverse(conj(outer))])
savefig("paper/figures/example_1_k23_srg_full_lfr.pdf") 

# compute SRGs for each individual operator
alphas = collect(-5:.01:5)
frequencies = exp10.(collect(-3:0.01:3))
phis = collect(0:0.001:pi)
plot()
plot_lti_srg!(Gyu, alphas, phis, frequencies)
savefig("paper/figures/example_1_k23_srg_Gyu.pdf") 

plot()
plot_lti_srg!(Gzu, alphas, phis, frequencies)
savefig("paper/figures/example_1_k23_srg_Gzu.pdf") 

plot()
plot_lti_srg!(Gyw, alphas, phis, frequencies)
savefig("paper/figures/example_1_k23_srg_Gyw.pdf") 

# compute rough gain bound using max gain of each operator, quicker but more conservative
maximum(abs.(srg_Gyw))*maximum(abs.(srg_return_upperhalf))*maximum(abs.(srg_Gzu))+maximum(abs.(srg_Gyu))
