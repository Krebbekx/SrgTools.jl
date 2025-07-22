using SrgTools, Plots, LaTeXStrings

# Nonlinear network: controlled Lur'e system with load side nonlinearity
# This example uses Julia functions instead of TransferFunction objects

# Define the stable controller and unstable LTI plant
K(s) = 1/(s+1)
G(s) = 3/((s-2)*(s/10+1))

# by Nyquist criterion, k1 > 2/3 is required
k1, k2 = .5, 1.5 # loop transformation gains. Make sure that these result in stable plants!

Gt(s) = G(s)/(1+k2*G(s)) # loop transformed G

S(s) = 1/(1+k1*K(s)*Gt(s)) # sensitivity function

# elements of the LFR matrix
G_z1u(s) = S(s)*K(s)
G_yu(s)  = S(s)*k1*K(s)*Gt(s)
G_z1w(s) = -S(s)*K(s)*Gt(s)
G_yw(s)  = S(s)*Gt(s)

Glfr(s) = [G_z1w(s) G_z1w(s) G_z1u(s); 
           G_yw(s)  G_yw(s)  G_yu(s) ;
           G_yw(s)  G_yw(s)  G_yu(s)]

Gzw(s) = Glfr(s)[1:2,1:2]
Gzu(s) = Glfr(s)[1:2,3:3]
Gyw(s) = Glfr(s)[3:3,1:2]
Gyu(s) = Glfr(s)[3:3,3:3]


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

# trick to plot the complex plane except for the inverted circle as the SRG of the inverted circle
farfield = 5
circle_inv_plot = [circle_inv; farfield; farfield+farfield*1im; -farfield+farfield*1im; -farfield-farfield*1im; farfield-farfield*1im; farfield]

plt = plot()
plot_lti_srg!(Gzw, alphas, phis, frequencies)
plot_srg!(circle_inv_plot)
xlims!((-3,3))
ylims!((-2.5,2.5))
annotate!(-1.4, -2, text(L"$\operatorname{SRG}(\Phi^\mathrm{a})^{-1}$", 26))
annotate!(.4, -.8, text(L"$\operatorname{SRG}(G_\mathrm{zw}^\mathrm{a})$", 26))
savefig("paper/figures/example_1_k0515.pdf") 
display(plt)

# compute the srg of (phi^{-1} - Gzw) in case of inner and outer radius
srg_Gzw_upperhalf = compute_lti_srg_boundary_upperhalfplane(Gzw, alphas, phis, frequencies)
circle_inv_upperhalf = [z for z in circle_inv if imag(z) >= 0]

# note that full srg, not only upperhalfplane, is needed for sum operation
srg_Gzw_chord = add_chords(srg_Gzw_upperhalf, 10) # add chords
srg_denom = [z1+z2 for z1 in circle_inv for z2 in -srg_Gzw_chord] # minus sign only with real numbers

# compute inner/outer radius of upper half plane part of denom. Only use inner radius.
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
savefig("paper/figures/example_1_k0515_srg_full_lfr.pdf") 

# compute rough gain bound
srg_Gzu = compute_lti_srg_boundary_upperhalfplane(Gzu, alphas, phis, frequencies)
srg_Gyw = compute_lti_srg_boundary_upperhalfplane(Gyw, alphas, phis, frequencies)
srg_Gyu = compute_lti_srg_boundary_upperhalfplane(Gyu, alphas, phis, frequencies)

maximum(abs.(srg_Gyw))*maximum(abs.(srg_return_upperhalf))*maximum(abs.(srg_Gzu))+maximum(abs.(srg_Gyu))

# case of SRG Gzw contained in Phi inverse

separation = minimum([minimum(abs.(srg_Gzw_upperhalf .- z)) for z in circle_inv])
# note that for stability, it is only required that the SRG is contained INSIDE a circle, defined by the INVERSE of the static NL.

# compute rough gain bound using max gain of each operator, quicker but more conservative
maximum(abs.(srg_Gyw))*(1/separation)*maximum(abs.(srg_Gzu))+maximum(abs.(srg_Gyu))