#Runge Kutta Stability functions
rksf(order)=rksffunction(z)=begin
    arr=[z^i/factorial(i) for i=0:order]
    return sum(arr)
end

g(order)=gn(r, i)=abs(rksf(order)(r.+i.*im))

xrange = -5:0.01:5
yrange = xrange

colors = ["magenta", "orange", "green", "red", "blue"]

p = contour(xrange, yrange, g(1), levels=[1], colorbar=false, aspect_ratio=1, color=colors[1])
for i in 2:5
    contour!(xrange, yrange, g(i), levels=[1], colorbar=false, aspect_ratio=1, color=colors[i])
end

savefig("01a.pdf")

#b)
#Verfahren 1
sf1(z)=1.0 + z + 1.0/2.0*z^2 + 1.0/6.0*z^3 + 1.0/12.0*z^4
#Verfahren 2
sf2(z)=1.0 + z + 1.0/2.0*z^2 + 1.0/6.0*z^3

gsf1(r, i)= abs(sf1(r+i*im))
gsf2(r, i)= abs(sf2(r+i*im))
p = contour(xrange, yrange, gsf1, levels=[1], colorbar=false, aspect_ratio=1, color=colors[1])
contour!(xrange, yrange, gsf2, levels=[1], colorbar=false, aspect_ratio=1, color=colors[2])
savefig("01b.pdf")

#d)
r_euler(z)=1.0 + z/(1.0-z)
r_middlepoint(z) = 1.0 + z/(1.0-(1.0/2.0)*z)
γ1 = (3.0+sqrt(3))/6
γ2 = (3.0-sqrt(3))/6
r_sdirk31(z) = 1.0 + (z + 1.0/2.0*z^2.0 - 2.0*z^2.0*γ1)/(1.0-z*γ1)^2.0
r_sdirk32(z) = 1.0 + (z + 1.0/2.0*z^2.0 - 2.0*z^2.0*γ2)/(1.0-z*γ2)^2.0

g_euler(r, i) = abs(r_euler(r+i*im))
g_middlepoint(r, i) = abs(r_middlepoint(r+i*im))
g_sdirk31(r, i) = abs(r_sdirk31(r+i*im))
g_sdirk32(r, i) = abs(r_sdirk32(r+i*im))

contour(xrange, yrange, g_euler, levels=[0,1], colorbar=false, aspect_ratio=1, color=colors[1])
savefig("01deuler.pdf")
contour(xrange, yrange, g_middlepoint, levels=[0,1], colorbar=false, aspect_ratio=1, color=colors[2])
savefig("01dmidpoint.pdf")
contour(xrange, yrange, g_sdirk31, levels=[0,1], colorbar=false, aspect_ratio=1, color=colors[3])
savefig("01dsdirk31.pdf")
contour(xrange, yrange, g_sdirk32, levels=[0,1], colorbar=false, aspect_ratio=1, color=colors[4])
savefig("01dsdirk32.pdf")
