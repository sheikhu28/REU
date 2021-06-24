include("ImaginaryTime.jl");

c = 1.1 * c₁;

u = u_ic[:,end];

function V(x,t)
    λ = 2;
    ϵ = 0.001;
    return λ*(1/sqrt(2pi*ϵ))*exp.(-((x .- c*t).^2)/(2ϵ));
end

n = 1; it = 1;
u_save = u;
tim = [0.0];
sd = s -> (im) .* ((vcat([s[2]+2*s[1]], diff(diff(s)), [s[Int(N-1)]-2*s[Int(N)]])/∂x^2) .- (V(x1,tim[end]+∂t) .* s) .- (abs.(s).^2 .* s) .+ s .- (G .* s));
while((n-1)*∂t<10)
    k1 = ∂t.*sd(u);
    k2 = ∂t.*sd(u.+(k1./2));
    k3 = ∂t.*sd(u.+(k2./2));
    k4 = ∂t.*sd(u.+(k3));
    u1 = u .+ ((k1 .+ (2 .* k2) .+ (2 .* k3) .+ k4)./6);
    global n += 1;
    if(mod(n,round(1/∂t))==0)
        global u_save = [u_save u1];
        push!(tim,(n-1)*∂t);
        global it += 1;
        # Plots.display(plot(tim,x1,abs.(u_save).^2));
        # sleep(0.1)
    end
    global u = u1; end


# Plots.display(heatmap(tim,x1,abs.(u_save).^2, c=:roma, camera=(60,45),ylim=(0,200)));
# png("1");

save("vals.jld","u",u_save,"x",x1,"t",tim);


# for i = 1 : size(u_save)[1]
#     Plots.display(plot(x1, abs.(u_save[:,i])));
# end
