using Plots
plotly()
include("../integration_methods.jl")

f(x, y_x) = -1000 * y_x + 1000 * sin(x) + cos(x)
df(x, y_x) = -1000
y(x)= ℯ.^(-1000.0 * x) + sin.(x)

hs = [1, 0.1, 0.01, 0.001, 0.00199, 0.002, 0.0021]
x = [collect(0:h:10) for h in hs]

y_comp = [y(x[i]) for i in 1:size(x,1)]

y_explicit = [explicit_euler_vec(f,1,0, hs[j], 10)  for j in 1:size(x,1)]

y_implicit = [implicit_euler_vec(f, df, 1, 0, hs[i], 10) for i in 1:size(x,1)]

plots = [plot(x[i],[y_comp[i], y_explicit[i], y_implicit[i]], title="h="*string(hs[i]), ylims=(-10,10)) for i in 1:size(x,1)]

plot(plots..., size=(1720,924))

savefig("Aufgabe3.html")

print("\n")
#
# y_explicit = [explicit_euler(f,1,0,1,i) for i in x]
#
# print(y_explicit)
