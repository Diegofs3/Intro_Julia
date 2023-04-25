#### Paquetes ####
#using Pkg
#Pkg.add("PyPlot")
#Pkg.add("FFTW")
#Pkg.add("LinearAlgebra")
using PyPlot,FFTW,LinearAlgebra

function wavevector(L,Nx);
    knx = collect(1:Nx) .-(Nx/2 +1)
    kx = 2pi*fftshift(knx)/L
end

function rk_s(tmax, dt, tplot, u0, kx, rko, pexp)
    Nx = length(kx)
    k1 = 1im*kx
    Lk = -k1.^3
    kmax = maximum(kx)
    filterk = zeros(Nx)
    filterk = exp.(-36.84*(abs.(kx)/kmax).^pexp)
    
    nplots = round(Int,tmax/tplot)
    plotgap = round(Int,tplot/dt)
    uarr = zeros(Nx,nplots)*1im
    
    u0 = u0*(1 + 0im) # Force complex initial function
    uhat0 = similar(u0)
    uhat = similar(u0)
    uhat1 = similar(u0)
    uhat2 = similar(u0)
    Px = plan_fft(u0) # compute FFTW plans
    
    for i = 1:nplots
        for n = 1:plotgap
            mul!(uhat0, Px, u0)
            for j = rko:-1:1
                uc = u0.*abs2.(u0)
                uhat1 = (uhat0 .+ dt*(Lk.*mul!(uhat, Px, u0) .+ 2im*mul!(uhat2, Px, uc))/j).*filterk
                ldiv!(u0, Px, uhat1)
            end
        end
        uarr[:, i] = u0
    end

    return uarr
end

Nx = 128
L = 8pi
h = L/Nx
x = (1:Nx)*h .- L/2
kx = wavevector(L, Nx)

dt  = 2.6/2000
A0 = 1e-2
uinit = 1 .+ A0*exp.(3im*x/4);

tmax = 120        # Final time 
rkor = 4        # Oder of the Runge-Kutta integrator
pexp = 36       # value for filter exponent. Use between 19 and 36 for best results
tplot = 0.1

uarr = rk_s(tmax, dt, tplot, uinit, kx, rkor, pexp);

############### Plots ##############

imshow(abs2.(uarr'), aspect="auto", origin="lower",
    extent=[-L/2, L/2, 0, tmax], interpolation = "bilinear", cmap=ColorMap("prism"), vmin=0, vmax=1)
ylabel(L"$Time$")
xlabel(L"$x$")
savefig("Graf8.svg")