include("../integration_methods.jl")


h = 0.1
t0 = 0
t_end = 5

x(t) = 1 ./144 .* t .^4.0 .+ π/2
x0 = [x(t0); sqrt(4/3 * (x(t0)-π/2)^(3/2))]


f(t,x) = [x[2]; sqrt(max(0.0, x[1]-π/2))]

t = collect(t0:h:t_end)

y = x(t)
y_euler = explicit_euler_vec(f, x0, t0, h, t_end)
y_rk = ERK_vec(f, x0, t0, h, t_end)

t0 = 0.01
x0 = [x(t0); sqrt(4/3 * (x(t0)-π/2)^(3/2))]
t_alt = t0:h:t_end

y_euler_alt = explicit_euler_vec(f, x0, t0, h, t_end)
y_rk_alt = ERK_vec(f, x0, t0, h, t_end)


using Plots
using LaTeXStrings
pyplot()

plot(t, y, label="Exact Solution", title=L"$\frac{1}{144}t^{4}+\frac{\pi}{2}$")
plot!(t, y_euler[:,1], label="Explicit Euler h="*L"$1\dot{}10^{-1}$")
plot!(t, y_rk[:,1], label="Runge Kutta h="*L"$1\dot{}10^{-1}$")

plot!(t_alt, y_euler_alt[:,1], label="Explicit Euler h="*L"$1\dot{}10^{-1}$"*" Verschoben um t="*L"$10^{-2}$")
plot!(t_alt, y_rk_alt[:,1], label="Runge Kutta h="*L"$1\dot{}10^{-1}$"*" Verschoben um t="*L"$10^{-2}$")

savefig("RK4_VS_ExplicitEuler.pdf")
