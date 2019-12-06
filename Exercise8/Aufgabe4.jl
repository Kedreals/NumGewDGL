include("../integration_methods.jl")
using LinearAlgebra


γ=6.672*10.0^-11;
mE=5.98*10.0^24;
mM=6.42*10.0^23;
mS=1.99*10.0^30;

# Gravity acting on the mars to the earth
FME(r)=(γ*mE*mM)*(r[1:2]-r[5:6])/(norm(r[1:2]-r[5:6])^3.0)
# Gravity acting on the earth to the sun
FES(r)=(γ*mE*mS)*(r[9:10]-r[1:2])/(norm(r[9:10]-r[1:2])^3.0)
# Gravity acting on the earth to the sun
FMS(r)=(γ*mM*mS)*(r[9:10]-r[5:6])/(norm(r[9:10]-r[5:6])^3.0)
f(t, r)= [  r[3:4];
            (-FME(r)+FES(r))./mE;
            r[7:8];
            (FME(r)+FMS(r))./mM;
            r[11:12];
            (-FES(r)-FMS(r))./mS];

rE0 = [150.0*10.0^9;0.0];
rM0 = [228.0*10.0^9;0.0];
rS0 = [0.0;0.0];
rdE0= [0.0; 29.0*10.0^3];
rdM0= [0.0; 24.0*10.0^3];
rdS0= [0.0; 0.0];

r0 = [rE0;rdE0;rM0;rdM0;rS0;rdS0];

t_erde=31557600;
t_mars=59355072;
h = 500;
t_end = t_mars

r_euler = explicit_euler_vec(f, r0, 0, h, t_end);
r_erk = ERK_vec(f, r0, 0,h, t_end);
r_ab5 = adams_bashforth_five(f,r0,0,h,t_end);

using Plots
pyplot()

# anim = @animate for i=2:100
#     plot(r_euler[(i-2)*trunc(Int,end/102)+1:i*trunc(Int, end/102),1], r_euler[(i-2)*trunc(Int,end/102)+1:i*trunc(Int, end/102),2], xlims=(-3*10.0^11,3*10.0^11), ylims=(-3*10.0^11,3*10.0^11), label="Euler Earth")
#     plot!(r_euler[(i-2)*trunc(Int,end/102)+1:i*trunc(Int, end/102),5], r_euler[(i-2)*trunc(Int,end/102)+1:i*trunc(Int, end/102),6], label="Euler Mars")
#     plot!(r_euler[(i-2)*trunc(Int,end/102)+1:i*trunc(Int, end/102),9], r_euler[(i-2)*trunc(Int,end/102)+1:i*trunc(Int, end/102),10], label="Euler Sun")
# end
# gif(anim, "Euler_Solar_System.gif", fps=15)

plot(r_euler[:,1], r_euler[:,2], label="Euler Earth")
plot!(r_erk[:,1],  r_erk[:,2],   label="RK Earth")
plot!(r_ab5[:,1],  r_ab5[:,2],   label="Adams Bashforth Earth")
plot!(r_euler[:,5], r_euler[:,6], label="Euler Mars")
plot!(r_erk[:,5],  r_erk[:,6],   label="RK Mars")
plot!(r_ab5[:,5],  r_ab5[:,6],   label="Adams Bashforth Mars")
plot!(r_euler[:,9], r_euler[:,10], label="Euler Sun")
plot!(r_erk[:,9],  r_erk[:,10],   label="RK Sun")
plot!(r_ab5[:,9],  r_ab5[:,10],   label="Adams Bashforth Sun")

savefig("Aufgabe4.pdf")
