include("../integration_methods.jl")
using LinearAlgebra

f(t, x)=(-t^2.0).*x
df(t, x) = -t^2.0 .* ones(size(x))
x0 = [ℯ]
t0 = 0
te = 1
hs = [2^float(-i) for i in 2:10]

x(t) = ℯ.^(1.0 .- (1.0/3.0).*(t.^3.0))

errors=zeros((16,size(hs,1)))

for (i,h) in enumerate(hs)
    t = collect(t0:h:te)
    x_exact = x(t)
    for j in 1:4
        if j == 1
            xstart = ERK_vec(f, x0, t0, h, t0+2*h, c=[0],A=[0],b=[1])
        elseif j == 2
            xstart = ERK_vec(f, x0, t0, h, t0+2*h, c=[0;1/2],A=[0 0;1/2 0],b=[0 1])
        elseif j == 3
            xstart = ERK_vec(f, x0, t0, h, t0+2*h, c=[0;1], A=[0 0; 1 0], b=[1/2 1/2])
        elseif j == 4
            xstart = ERK_vec(f, x0, t0, h, t0+2*h)
        end
        xstart = xstart'

        x_adams_bashford = adams_bashforth_three(f, xstart, t0, h, te)
        x_adams_moulton = adams_moulton_three(f, df, xstart, t0, h, te)
        x_nyström = Nyström_three(f, xstart, t0, h, te)
        x_milne_simpson = milne_simpson_two(f, df, xstart, t0, h, te)

        errors[(j-1)*4+1,i] = x_adams_bashford[end]-x_exact[end]
        errors[(j-1)*4+2,i] = x_adams_moulton[end]-x_exact[end]
        errors[(j-1)*4+3,i] = x_nyström[end] - x_exact[end]
        errors[(j-1)*4+4,i] = x_milne_simpson[end] - x_exact[end]
    end
end

plots=[]
for i in 1:4
    if i == 1
        titel = "Adams Bashford"
    elseif i == 2
        titel = "Adams Moulton"
    elseif i == 3
        titel = "Nyström"
    else
        titel = "Milne Simpson"
    end
    plt = plot(1:size(hs,1)-1, ooc(hs, errors[i,:]), label="explicit euler", title = titel, xlabel = "index of hs")
    plot!(1:size(hs,1)-1, ooc(hs, errors[4+i,:]), label="Runge Kutta 2")
    plot!(1:size(hs,1)-1, ooc(hs, errors[8+i,:]), label="Heun")
    plot!(1:size(hs,1)-1, ooc(hs, errors[12+i,:]), label="Runge Kutta 4")

    push!(plots, plt)
end

plot(plots..., size=(1920, 1040))
savefig("Aufgabe1.pdf")
