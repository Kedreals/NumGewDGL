include("../integration_methods.jl")

h=0.01
x0=[1;2]

#1.)
f(A)=(t,x)->A*x

X1(t)= ((1/2 + im).*ℯ.^((-1-5*im) .*t))' .*[1;im] + ((1/2- im).*ℯ.^((-1+5*im) .* t))' .*[1;-im]
X2(t)= ((1/2 + 1/2*im).*ℯ.^(-2 .*im.*t))' .*[1;-2*im] + ((1/2 - 1/2*im).*ℯ.^(2 .*im.*t))' .*[1;(2*im)]

c1 = [0;1/2;1/2;1]
c2 = [0;1/2;1;1]

A1 = [0   0   0 0;
      1/2 0   0 0;
      0   1/2 0 0;
      0   0   1 0];

A2 = [0   0   0 0;
      1/2 0   0 0;
      1/4 1/4 0 0;
      0   -1  2 0];

b1 = [1/6 1/3 1/3 1/6];
b2 = [1/6 0 2/3 1/6];

t0=0
t_end=4
f1 = f([-1 5; -5 -1])
exact1 = X1(collect(t0:h:t_end))
euler1 = explicit_euler_vec(f1, x0, t0, h, t_end)
rk41 = ERK_vec(f1, x0, t0, h, t_end, c=c1,A=A1, b=b1)
england1 = ERK_vec(f1, x0, t0, h, t_end, c=c1,A=A2, b=b2)
r381 = ERK_vec(f1, x0, t0, h, t_end, c=c1,A=A2, b=b2)
rk41e = rk41[end,:]-exact1[:,end]
england1e = england1[end,:]-exact1[:,end]
r381e = r381[end,:]-exact1[:,end]

print(
"errorsquares of first problem
    RK4:     "*string(real(sum(rk41e.^2)))*"
    England: "*string(real(sum(england1e.^2)))*"
    3/8 Rule:"*string(real(sum(r381e.^2)))*"\n")

t0=0
t_end=15
f2 = f([0 1; -4 0])
euler2 = explicit_euler_vec(f2, x0, t0, h, t_end)
h = 0.1
exact2 = X2(t0:h:t_end)
rk42 = ERK_vec(f2, x0, t0, h, t_end, c=c1,A=A1, b=b1)
england2 = ERK_vec(f2, x0, t0, h, t_end, c=c1,A=A2, b=b2)
r382 = ERK_vec(f2, x0, t0, h, t_end, c=c1,A=A2, b=b2)

rk42e = rk42[end,:]-exact2[:,end]
england2e = england2[end,:]-exact2[:,end]
r382e = r382[end,:]-exact2[:,end]

print(
"errorsquares of first problem
    RK4:     "*string(real(sum(rk42e.^2)))*"
    England: "*string(real(sum(england2e.^2)))*"
    3/8 Rule:"*string(real(sum(r382e.^2)))*"\n")

using Plots
using LaTeXStrings
pyplot()

subplots1=[]
p = plot(real.(exact1[1,:]), real.(exact1[2,:]), label="Exact Solution")
plot!(euler1[:,1], euler1[:,2], label="Euler")
push!(subplots1, p)
p = plot(real.(exact2[1,:]), real.(exact2[2,:]), label="Exact Solution")
plot!(euler2[:,1], euler2[:,2], label="Euler")
push!(subplots1, p)
plot(subplots1..., size=(1024, 768))
savefig("Aufgabe1.pdf")


subplots2=[]
p = plot(real.(exact1[1,:]), real.(exact1[2,:]), label="Exact Solution")
plot!(rk41[:,1], rk41[:,2], label="RK4")
plot!(england1[:,1], england1[:,2], label="England Regel")
plot!(r381[:,1], r381[:,2], label=L"$\frac{3}{8}$"*" Rule")
push!(subplots2, p)
p = plot(real.(exact2[1,:]), real.(exact2[2,:]), label="Exact Solution")
plot!(rk42[:,1], rk42[:,2], label="RK4")
plot!(england2[:,1], england2[:,2], label="England Regel")
plot!(r382[:,1], r382[:,2], label=L"$\frac{3}{8}$"*" Rule")
push!(subplots2, p)
plot(subplots2..., size=(1200, 900))
savefig("Aufgabe2.pdf")
