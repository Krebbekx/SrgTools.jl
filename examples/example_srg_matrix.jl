using SrgTools, Plots

# This script shows how to the LTI SRG algorithm can be used to compute the SRG of a matrix. In this case, the matrix is simply modeled as a constant transfer function. Therefore, we can take a frequency vector of only one value. 

frequencies = [1] # only one value (it does not matter what value) 
alphas = collect(-5:.1:5)
phis = collect(0:0.01:pi)

Gzw(s) = [0 5 0; 0 0 -1; 1 1 0]; # constant matrix function example
#Gzw(s) = [0 0 0; 0 0 -1; 1 1 0]; # constant matrix from FF setup

# animation of the max (blue) and min (blue) gain circles
# plt = plot()
# plot_lti_srg_circles!(Gzw, alphas, phis, frequencies)

# plot of the resulting SRG
plt = plot()
plot_lti_srg!(Gzw, alphas, phis, frequencies, :blue, 0.7)