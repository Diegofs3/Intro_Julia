#### Paquetes ####
#using Pkg
#Pkg.add("PyPlot")
#Pkg.add("FFTW")
using PyPlot,FFTW

####### Rk4 Improved ######

function rk104(u,dt,g,ep)
    q1 = u 
    q2 = u
    for i =1:5
        #q1 = q1 .+ dt* Frh(q1,g,ep)/6
        #Linea equivalente
        q1 .+= dt* Frh(q1,g,ep)/6
    end
    q2 = (q2 .+ 9*q1)/25
    q1 = 15*q2 .- 5*q1
    for i =6:9
        #q1 = q1 .+ dt* Frh(q1,g,ep)/6
        q1 .+= dt* Frh(q1,g,ep)/6
    end
    q1 = q2 .+ 3*q1/5 .+ dt*Frh(q1, g, ep)/10
    return u = q1 
end

##### Aplicado en el problema de Korteweg #####

function Frh(u, g, ep)
    qrh = ifft(ep.*u)
    qrh = g.*fft(qrh.^2)
end

function wavevector(L,Nx);
    knx = collect(1:Nx) .-(Nx/2 +1)
    kx = 2pi*fftshift(knx)/L
end

N = 256                           #Numero de celdas del dominio, recomendable 2^n
L = 4pi                             #Longitud del analisis
h = L/N                            #Ancho de cada celda
dt = 0.04/N^2
x = (0:N-1)*h .-L/2
k = wavevector(L,N)
kmax = maximum(k)
fillterk = exp.(-36.84*(abs.(k/kmax)).^36)
A , B = 25,16
u0 = 3*A^2*sech.(A*(x .+ 2)/2).^2 .+ 3*B^2*sech.((B*(x .+1))/2).^2
tf = 0.006
nt = Int(div(tf,dt))+1
uarr = zeros(N,nt+1)
uarr[:,1] = u0
ik3 = 1im*dt*k.^3
ihk = -0.5im*k
ep = exp.(ik3);
nj = 1
Uhat = fft(u0)
for t = 0:dt:tf
    global nj , Uhat
    global ep = exp.(ik3* nj)
    g = ihk./ ep
    Uhat = rk104(Uhat,dt,g,ep).*fillterk
    nj += 1
    uarr[:,nj] = real(ifft(ep.*Uhat))
end
for i in [60,300,450]
    plot(uarr[:,i])
end
savefig("Graf7.svg")
clf()
imshow(uarr',aspect="auto",origin="lower",extent=[-L/2-4,L/2+4,-0.0001,tf], interpolation = "bicubic",cmap = ColorMap("jet"))
colorbar()
savefig("Graf8.svg")

