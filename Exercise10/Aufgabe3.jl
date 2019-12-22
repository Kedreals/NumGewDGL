include("../integration_methods.jl")

f(t, x) = -x .+ ℯ^(-t)*cos(t)
df(t,x) = -1.0

x(t) = ℯ.^(-t).*sin.(t)

x0 = 0
t0 = 0
h = 10.0^(-3.0)
t_end = 10

αs = [-1,-0.99999,-0.9, 0.9, 0.99999, 1]

x0s = ERK_vec(f, [x0], t0, h, t0+h)

ts = collect(t0:h:t_end)

x_exact = x(ts)
pl = plot(ts, x_exact, label="Analytical Solution", linestyle= :dot, linewidth= 2)
for α in αs
    x_alpha = special_linear_2_schritt_verfahren(f, df, x0s, t0, h, t_end, α, newton_iterationen=5)
    plot!(ts, x_alpha, label="Alpha=$α",linestyle= :dash)
end

display(pl)
