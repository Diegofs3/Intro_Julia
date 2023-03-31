######## Derivadas? FFT ##########
using Pkg
#Pkg.add("PyPlot")
#Pkg.add("FFTW")
using PyPlot, FFTW

######### Funcion que no se porque está realmente #######

function wavevector(L,Nx);
    knx = collect(1:Nx) .-(Nx/2 +1)
    kx = 2pi*fftshift(knx)/L
end

######### Condiciones Iniciales ##########

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

plot(x,f)
plot(x,fs)
savefig("Graf6.svg")
clf()
print("Eso es todo \n")