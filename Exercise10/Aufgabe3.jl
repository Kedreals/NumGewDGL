include("../integration_methods.jl")

f(t, x) = -x .+ ℯ^(-t)*cos(t)
df(t,x) = -1.0

x(t) = ℯ.^(-t).*sin.(t)

x0 = 0
t0 = 0
h = 5*10.0^(-1.0)
t_end = 10

αs = [-1,-0.99999,-0.9, 0.9, 0.99999, 1]

x0s = ERK_vec(f, [x0], t0, h, t0+h)

ts = collect(t0:h:t_end)

pl = plot(t0:0.001:t_end, x(t0:0.001:t_end), label="Analytical Solution", linestyle= :dot, linewidth= 2)
for α in αs
    x_alpha = special_linear_2_schritt_verfahren(f, df, x0s, t0, h, t_end, α, newton_iterationen=5)
    plot!(ts, x_alpha, label="Alpha=$α",linestyle= :dash)
end

savefig("Aufgabe3.pdf")
