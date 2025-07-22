using SrgTools, Plots

G(s::Number) = [s/(s+1) 2*s/(s+1); s^2/(s+1)^2 s/(s+1)^2] # 2x2 transfer function

frequencies = exp10.(collect(-3:0.01:3)) # grid of frequency points to evaluate over

alphas = collect(-5:.05:5) # real axis points to compute the SRG radius with 

phis = collect(0:0.01:pi) # angular resolution of the outer bound of the SRG in the upper half plane

# to obtain the complex points that parameterize the boundary (as a polygon) of the SRG in the upper half complex plane
srg_upperhalf = compute_lti_srg_boundary_upperhalfplane(G, alphas, phis, frequencies)

plt = plot()
plot_srg!(srg_upperhalf)
xlims!((-3,3))
ylims!((-3,3))
display(plt)

srg_upperhalf = [srg_upperhalf; reverse(conj(srg_upperhalf))]

srg2 = [5+1*exp(1im*phi)  for phi in collect(range(-pi, stop=pi, length=100))] # circle arc SRG

srg_sum = [z1+z2 for z1 in srg_upperhalf for z2 in srg2]

inner_rad, outer_rad = inner_outer_radius_srg_binning(srg_sum .- 3, 50);

srg = [inner_rad .+ 3; reverse(outer_rad .+ 3)]

plt = plot()
plot!(real(srg), imag(srg))
plot!(real(srg_sum), imag(srg_sum))
display(plt)