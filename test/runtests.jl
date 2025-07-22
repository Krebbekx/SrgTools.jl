using Test, SrgTools

########## Stuff to make into tests later ##########

# test it, seems to work!
# G22(s::Number) = [s/(s+1) s/(s+1); s/(s+1)^2 s/(s+1)^2] # 2x2 TF

# G12(s::Number) = [s/(s+1) s/(s+1)] # 1x2 TF

# G21(s::Number) = [s/(s+1); s/(s+1)^2 ] # 2x1 TF

# G11(s::Number) = [1/(s^2+s+1)] # 1x1 TF

# Gs(s::Number) = 1/(s^2+s+1) # scalar TF

# frequencies = [1, 2, 3]

# alphas = [-5, -1, 0, 2, 10]

# compute_srg_circles_lti(G22, alphas, frequencies)

# compute_srg_circles_lti(G12, alphas, frequencies)

# compute_srg_circles_lti(G21, alphas, frequencies)

# compute_srg_circles_lti(G11, alphas, frequencies)

# compute_srg_circles_lti(Gs, alphas, frequencies)


# frequencies = exp10.(collect(-3:0.01:3))

# alphas = collect(-5:.05:5)

# phis = collect(0:0.01:pi)

# plot()
# plot_lti_srg!(G22, alphas[1:end], phis, frequencies)

# xlims!((-2,2))
# ylims!((-2,2))

# plot_lti_srg_old!(Gs, alphas[1:end], phis, frequencies) # compare/check
# plot!(xaxis=true)
# savefig("myplot.pdf")
# srg_upperhalf = compute_lti_srg_boundary_upperhalfplane(G22, alphas, phis, frequencies)
# srg = [srg_upperhalf; conj(reverse(srg_upperhalf))]
# boundary = srg_plus_circle_boundary(srg, 2, 0, 100)

# plot(real(boundary), imag(boundary))


########## testing stuff ##########


# plot()
# phis = collect(0:0.02:2.01*pi)
# plot_lti_srg_circles!(phis, alphas, max_radii)
# plot_lti_srg_circles!(phis, alphas, min_radii)


# function test_varphi_finder(alpha, radius, z0, phis)
#     varphis = [compute_varphi_srg_center(z0, alpha, radius, phi) for phi in phis]

#     z1, z2 = zeros(Complex, length(phis)), zeros(Complex, length(phis))
#     for i = 1:length(phis)
#         z1[i] = alpha + radius*exp(1im*varphis[i])
#         z2[i] = z0 + abs(alpha + radius*exp(1im*varphis[i])-z0)*exp(1im*phis[i])
#     end
#     # z1 = [alpha + radius*exp(1im*phi) for phi in phis]
#     # z2 = [z0 + abs(alpha + radius*exp(1im*varphi)-z0)*exp(1im*varphi) for varphi in varphis]

#     display(plot(real(z1),imag(z1)))
#     display(plot!(real(z2),imag(z2)))

#     print(maximum(abs.(z1-z2)))
# end 

# test_varphi_finder(-5, 6, -1, collect(0:0.01:pi))

# """
# Function to test the dimension cast in compute_srg_circles_lti
# """
# function test_dimensions_G(G)
#     if length(size(G(0))) == 2 # matrix and row vector case
#         q, p = size(G(0))
#     elseif length(size(G(0))) == 1 # column vector case, or scalar in []
#         p = 1
#         q = size(G(0))[1]
#     else # scalar case
#         p = q = 1
#     end

#     # correct for dimensions and create correct identity
#     G_alpha(s) = [G(s); zeros(max(p-q,0),p)]

#     return G_alpha(0)
# end