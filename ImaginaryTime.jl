using Plots, LinearAlgebra
gr()

s = 0.4;
∂x = 1/30;
∂t = s * ∂x^2;

ll = 300;
N = Int(round(2*ll/∂x) - 1);
x = range(-ll, ll, length = Int(N)+2);
x1 = zeros(length(x)-2);
for i = 2:N+1
    x1[Int(i-1)] = x[Int(i)];
end

x₀ = 0.69669;
c = 0.3;

g(x,a) = sqrt((c^2/2)+(1-(c^2/2))*(tanh(sqrt(0.5-(c^2)/4)*(x+a)))^2);
f(x) = atan(sqrt(2c^2-c^4)/(exp(sqrt(2-c^2)*(x+x₀))+c^2-1));

function R(x)
    if(x<0)
        return g(x,-x₀);
    else
        return g(x,x₀);
    end
end

function ϕ(x)
    if(x<0)
        return 2*f(0)-f(-x);
    else
        return f(x);
    end
end
w = 250;
γ(x) = (1-exp(-((x)/w)^12));
G = 10 .* γ.(x1);

u = (R.(x1) .* exp.(im .* ϕ.(x1)));

function v(x)
    λ = 2;
    ϵ = 0.0001;
    return λ*(1/sqrt(2pi*ϵ))*exp.(-(x.^2)/(2ϵ));
end
Vx = v(x1);

μ = 0;

∂t2 = -(im)*∂t; #imaginary time

n = 1; it = 1;
u_ic = u;
tim = [0.0];
err = [0.0];
sd = s ->  ((c .* vcat([s[2]],[s[i+1]-s[i-1] for i = 2:N-1], [s[N-1]]))/2∂x) .+ ((im) .* ((vcat([s[2] - 2*s[1]], diff(diff(s)), [s[Int(N-1)] - 2*s[Int(N)]])/2∂x^2) .- (Vx .* s) .- (abs.(s).^2 .* s) .+ s .- (G .* s) .+ (μ .* s))) ;
while(true)
    k1 = ∂t2.*sd(u);
    k2 = ∂t2.*sd(u.+(k1./2));
    k3 = ∂t2.*sd(u.+(k2./2));
    k4 = ∂t2.*sd(u.+(k3));
    u1 = u .+ ((k1 .+ (2 .* k2) .+ (2 .* k3) .+ k4)./6);
    global n += 1;
    if(mod(n,Int(1/∂t))==0)
        global u_ic = [u_ic u1];
        push!(tim,(n-1)*∂t);
        push!(err, norm(abs.(u1) .- abs.(u)));
        global it += 1;
        Plots.display(plot(tim,err,ylims=(10^-8,0.005),yscale=:log10));
        # sleep(0.1)
    end
    u2 = u;
    global u = u1;
    if(norm(abs.(u1).^2 .- abs.(u2).^2) < 10^-8)
        break;
    end end


Plots.display(plot(x1,abs.(u_ic[:,end]).^2));

# for i = 1 : size(u_save)[1]
#     Plots.display(plot(x1, abs.(u_save[:,i])));
# end
