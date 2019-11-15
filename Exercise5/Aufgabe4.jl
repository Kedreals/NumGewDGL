using Plots
using LaTeXStrings
pyplot()
include("../integration_methods.jl")

ρ = 998.2;
σ = 0.07274;
g = 9.81;
z0 = -1.0 * 10.0^-3;
kss = [0.2, 0.3, 0.4, 0.5] .* 10.0^3;
r0 = 0;
φ0 = 0;

F(ks) = (s, x) -> [cos.(x[3]); sin.(x[3]); ((x[1] > 0) ? 2.0*ks .- (ρ*g/σ).*x[2].-sin.(x[3])./x[1] : ks .- (ρ*g)/(2.0*σ) .* x[2])]

x_0 = [r0; z0; φ0];

h = 0.1 * 10^-3

vals = RK3_vec(F(kss[1]), x_0, 0, h, 0.005);
ind = vals[:,2] .<= 0
plot(vals[ind,1], vals[ind,2], label="ks="*string(kss[1])*L"$m^{-1}$", xlabel=L"m", ylabel=L"m", title=L"Runge Kutta 3")
for i in 2:size(kss,1)
    vals = RK3_vec(F(kss[i]), x_0, 0, h, 0.005)
    ind = vals[:,2] .<= 0
    plot!(vals[ind,1],vals[ind,2], label="ks="*string(kss[i])*L"$m^{-1}$")
end

savefig("Tropfen_Euler_RK3.pdf")
