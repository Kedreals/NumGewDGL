include("../integration_methods.jl")

f(t,x)=[x[2]; (-10.0*(x[1]^2.0 - 1.0)*x[2] - x[1])]
x0 = [1.0; 1.0]
t0 = 0

x_h, t_h = ARK_vec(f, x0, t0, 10^(-2), 10)

plot(t_h, [x_h[i][1] for i in 1:size(t_h,1)], xticks=1:10)
savefig("Van_Der_Pol.pdf")

hs = t_h[2:end]-t_h[1:end-1]

plot(t_h, [hs...; hs[end]], xticks=1:10)
savefig("hs.pdf")
