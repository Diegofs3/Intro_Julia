######## ¿Derivadas e integrales? FFT ##########
#using Pkg
#Pkg.add("PyPlot")
#Pkg.add("FFTW")
using PyPlot, FFTW

######### Funcion que no se porque está realmente #######

function wavevector(L,Nx);
    knx = collect(1:Nx) .-(Nx/2 +1)
    kx = 2pi*fftshift(knx)/L
end

########### Derivada?? ######

function ffteriv(f,k1)
    dfk = k1.*fft(f)
    df = real(ifft(dfk))

########## Integral de rk??? ########

function rk_s(tf,dt,u0,kx,rko)
    k1 = lim*kx
    u1 = u0
    for t = 0:dt:tf
        for j = rko:-1:1
            dudx = ffteriv(u1.k1)
            u1 = u0 .+ dt*cx.*dudx/j
        end
        u0 = real(u1)
    end
    return u0
end



######### Condiciones Iniciales derivada ##########

Nx = 2048                           #Numero de celdas del dominio, recomendable 2^n
L = 8pi                             #Longitud del analisis
h = L/Nx                            #Ancho de cada celda
x = (0:Nx-1)*h                      #Valores de x por celda
x0 = 10                             #Cte arbitraria
sig = 2                             #Ancho Gaussiana
kx = wavevector(L,Nx)               #Valores de k para F{f(x)} = g(k)
ct = 5                              #Desplazamiento de la Gausiana??          
f = exp.(-(x .-x0).^2/sig^2)        #Curva Gaussiana
fk = fft(f)                         #Transformada de una Gaussiana = Gaussiana
dev_f =real(ifft(im*kx.*fk))        #
z = exp.(-im*ct*kx)                 #Valores de g(k)
fs = real(ifft(z.*fk))              #Retorno al espacio de configuraciones

#Graficos

plot(x,f)
plot(x,fs)
savefig("Graf6.svg")
clf()

########## C.I Integral ############

Nx = 128
L0 = 2pi
h = L/Nx
x = (1:Nx)*h
kx = wavevector(L,Nx)
cx = 0.2 + sin(x-1)^2
dt = 0.001
u0 = exp.(-100*(x-1)^2)

###### Completar y ojala hacer de nuevo a gusto ######

print("Eso es todo \n")