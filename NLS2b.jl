# include("ImaginaryTime.jl")

u = u_ic[:,end];

γ(x,t) = γ(x+c*t);

n = 1; it = 1;
u_save = u;
tim = [0.0];
while((n-1)*∂t<20)
    G2 = γ.(x1,tim[end]+∂t);
    local sd = s ->  ((c .* vcat([s[2]],[s[i+1]-s[i-1] for i = 2:N-1], [-s[N-1]]))/2∂x) .+ ((im) .* ((vcat([s[2] - 2*s[1]], diff(diff(s)), [s[Int(N-1)] - 2*s[Int(N)]])/∂x^2) .- (Vx .* s) .- (abs.(s).^2 .* s) .+ s .- (G2 .* s)));
    k1 = ∂t.*sd(u);
    k2 = ∂t.*sd(u.+(k1./2));
    k3 = ∂t.*sd(u.+(k2./2));
    k4 = ∂t.*sd(u.+(k3));
    u1 = u .+ ((k1 .+ (2 .* k2) .+ (2 .* k3) .+ k4)./6);
    global n += 1;
    if(mod(n,25)==0)
        global u_save = [u_save u1];
        push!(tim,(n-1)*∂t);
        global it += 1;
    end
    global u = u1; end

Plots.display(heatmap(tim,x1,abs.(u_save).^2, c=:roma,camera=(90,45),ylims=(-50,50)));

# u_shifted = u_save[:,1];
# for i = 2:size(u_save)[2]
#     u_temp = copy(u_save[:,1]);
#     step = Int(round(i*c*25*∂t/∂x));
#     j = 1+step;
#     while(j<=size(u_save)[1])
#         global u_temp[j] = u_save[j-step,i];
#         global j+=1;
#     end
#     global u_shifted = [u_shifted u_temp];end
#
# Plots.display(surface(tim,x1,abs.(save .- u_shifted) .^2,c=:roma));
# plot!(camera=(90,45));

# for i = 1 : size(u_save)[1]
#     Plots.display(plot(x1, abs.(u_save[:,i])));
# end
