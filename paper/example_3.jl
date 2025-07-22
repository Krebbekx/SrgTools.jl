# Example comparison with Timo 2025 MIMO paper

using SrgTools, Plots, LaTeXStrings, ControlSystemsBase
gr() # use GR backend for plotting 

s = tf("s") # define the Laplace variable

P11 = 0.1/(s+1)
P12 = 1/(s^3+5*s^2+2*s+1)
P21 = P12/10
P22 = 0.2/(s+5)

P = [P11 P12; P21 P22]

Glfr = [P11 P12 P11 P12;
        P21 P22 P21 P22;
        P11 P12 P11 P12;
        P21 P22 P21 P22]

@assert(maximum(real(pole(Glfr)))<0) # check stability

Gzw = Glfr[1:2,1:2]
Gzu = Glfr[1:2,3:4]
Gyw = Glfr[3:4,1:2]
Gyu = Glfr[3:4,3:4]


frequencies = exp10.(collect(-5:0.01:5)) # logspace frequencies
alphas = collect(-5:.01:5)
phis = collect(0:0.005:pi)

# build static NL circle
radius = sqrt(0.1)
angles_circle = collect(range(0, stop=2*pi, length=300))
radius_circle = sqrt(0.1)
center_circle = 0
circle = center_circle .+ radius_circle .* exp.(1im*angles_circle)
circle_inv = invert_srg(circle)

# trick to plot the complex plane except for the inverted circle as the SRG of the inverted circle
farfield = 5 # for SRG plot
circle_inv_plot = [circle_inv; farfield; farfield+farfield*1im; -farfield+farfield*1im; -farfield-farfield*1im; farfield-farfield*1im; farfield]


plt = plot()
plot_srg!(circle_inv_plot)
plot_lti_srg!(Gzw, alphas, phis, frequencies)
xlims!((-4,4))
ylims!((-4,4))
annotate!(-1.7, 3, text(L"$\operatorname{SG}_0(\Phi)^{-1}$", 22))
annotate!(.5, -0.9, text(L"$\operatorname{SG}_0(G)$", 22))
# savefig("scripts/figures/two_masses_Gzw_Phi.pdf") 
display(plt)


# look at the individual SRG Plots
plot()
plot_lti_srg!(Gyu, alphas, phis, frequencies)
# savefig("scripts/figures/example_2_srg_Gyu.pdf") 

plot()
plot_lti_srg!(Gzu, alphas, phis, frequencies)
# savefig("scripts/figures/example_2_srg_Gzu.pdf") 

plot()
plot_lti_srg!(Gyw, alphas, phis, frequencies)
# savefig("scripts/figures/example_2_srg_Gyw.pdf") 


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
maximum(abs.(srg_upperhalf))

plot()
plot_srg!([outer; reverse(conj(outer))])
# savefig("scripts/figures/example_2_srg_full_lfr.pdf") 
srg_bound1 = maximum(abs.(srg_upperhalf))


srg_bound2 = maximum(abs.(srg_Gyu)) + maximum(abs.(srg_Gzu))*maximum(abs.(srg_Gyw))*maximum(abs.(srg_return))  # worst case bound

srg_bound1 < srg_bound2


#### direct application of incremental homotopy theorem! Yields better bound


srg_P_upperhalf = compute_lti_srg_boundary_upperhalfplane(P, alphas, phis, frequencies)

srg_P_upperhalf_inv = invert_srg(srg_P_upperhalf)

circle_upperhalf = [z for z in circle if imag(z) >= 0]

plt = plot()
plot_lti_srg!(P, alphas, phis, frequencies)
xlims!((-1,1))
ylims!((-2,2))
annotate!(-0, 1.5, text(L"$\operatorname{SRG}(P)$", 22))
xticks!([-0.5,0,0.5])
savefig("paper/figures/example_3_P.pdf") 
display(plt)

plt = plot()
plot_srg!(circle)
plot_srg!(srg_P_upperhalf_inv)
plot_srg!(conj(srg_P_upperhalf_inv))
xlims!((-5,4))
ylims!((-3,3))
annotate!(-1.7, 0, text(L"$\operatorname{SG}_0(\Phi)$", 22))
annotate!(-1, 2, text(L"$\operatorname{SRG}(P \;)^{-1}$", 22))
savefig("paper/figures/example_3_denom.pdf") 
display(plt)

srg_denom = [z1+z2 for z1 in circle for z2 in invert_srg(srg_P_upperhalf)]

srg_bound = 1/minimum(abs.(srg_denom))

inner, outer = improved_chord_sum(circle_upperhalf, invert_srg(srg_P_upperhalf), N, n_bins);

srg_upperhalf = [invert_srg(inner); reverse(invert_srg(outer))];
maximum(abs.(srg_upperhalf))

plot()
plot_srg!([invert_srg(inner); reverse(conj(invert_srg(inner)))])
annotate!(-1, 2, text(L"$\operatorname{SG}_0(H_1)$", 22))
xlims!((-2,2))
savefig("paper/figures/example_3_H1.pdf") 
display(plt)

srg_bound1 = maximum(abs.(srg_upperhalf))


##### use incremental homotopy on full closed loop

H2 = [1.7/(s^2+2*s+1) 0; 0 1.7/(s^2+3*s+3)] 

G = tf(feedback(ss(P), ss(H2)));

plot()
plot_srg!(circle_inv_plot)
plot_lti_srg!(G, alphas, phis, frequencies)

maximum(real(pole(G)))<0


srg_G_upperhalf = compute_lti_srg_boundary_upperhalfplane(G, alphas, phis, frequencies)

circle_upperhalf = [z for z in circle if imag(z) >= 0]

srg_denom = [z1+z2 for z1 in circle for z2 in invert_srg(srg_G_upperhalf)]

srg_bound = 1/minimum(abs.(srg_denom))

inner, outer = improved_chord_sum(circle_upperhalf, invert_srg(srg_G_upperhalf), N, n_bins);

srg_upperhalf = [invert_srg(inner); reverse(invert_srg(outer))];
maximum(abs.(srg_upperhalf)) # gain bound of 2.07 is obtained

plot()
plot_srg!([invert_srg(inner); reverse(conj(invert_srg(inner)))])
xlims!((-2,2))
annotate!(-1.2, 1.2, text(L"$\operatorname{SG}_0(T)$", 22))
savefig("paper/figures/example_3_T.pdf") 
display(plt)

srg_bound1 = maximum(abs.(srg_upperhalf))