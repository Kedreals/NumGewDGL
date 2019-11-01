using Plots
pyplot()
include("../integration_methods.jl")

f(x, y_x) = -1000 * y_x + 1000 * sin(x) + cos(x)
df(x, y_x) = -1000
y(x)= ℯ.^(-1000.0 * x) + sin.(x)

hs = [1, 0.1, 0.01, 0.001, 0.00199, 0.002, 0.0021]
x = [collect(0:h:10) for h in hs]

y_comp = [y(x[i]) for i in 1:size(x,1)]

y_explicit = [explicit_euler_vec(f,1,0, hs[j], 10)  for j in 1:size(x,1)]

y_implicit = [implicit_euler_vec(f, df, 1, 0, hs[i], 10) for i in 1:size(x,1)]

plots = []
for i in 1:size(x,1)
    p = plot(x[i], y_comp[i], label="reference", title="h="*string(hs[i]), ylims=(-10,10))
    plot!(x[i], y_explicit[i], label="explicit")
    plot!(x[i], y_implicit[i], label="implicit")
    push!(plots, p)
end

plot(plots..., size=(1720,924))

savefig("Aufgabe3.pdf")
