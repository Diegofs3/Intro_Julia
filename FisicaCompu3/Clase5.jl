##### Pendulo Forzado?? ######
using Pkg
#Pkg.add("FFTW")
using PyPlot,FFTW

##w0 frecuencia motor normalizado, b intensidad fuerza impulsora, q roce, y angulo
function Pendulo(y,q,b,w0,t)    
    g = [y[2],-q*y[2]- sin(y[1]) + b*cos(w0*t)]
end
function RK4(y,t,dt,q,b,w0)     #RungeKutta4
    n = length(y)      
    hk1 = zeros(n)
    hk2 = similar(hk1)
    hk3 = similar(hk1)
    hk4 = similar(hk1)

    hk1 = Pendulo(y , q , b , w0 , t )
    hk2 = Pendulo(y .+ dt*hk1/2 , q , b , w0 , t + dt/2)
    hk3 = Pendulo(y .+ dt*hk2/2 , q , b , w0 , t + dt/2)
    hk4 = Pendulo(y .+ dt*hk3 , q , b , w0 , t + dt)
    
    resp = y .+ dt*(hk1 .+ 2*(hk2 .+ hk3) .+hk4)/6
    resp[1] = resp[1] - 2pi*(abs(resp[1]) > pi)*abs(resp[1])/resp[1]    #Esto ajusta un -pi a pi cuando el pendulo pasa del punto más alto
    return resp                                                         #Igual tengo que investigar como funciona ese condicional en el codigo
end

#Definimos valores
nt = 2^16               #Numero de pasos
w0= 2/3                 #frecuencia motor externo
q= 0.5                  #Roce/amortiguamiento
b= 0.9                  #Intensidad de la fuerza externa
fp = 192                #ni idea
dt = (2pi/w0)/fp        #Salto temporal
tfin = (nt-1)*dt        #Ultima iteración
y = zeros(2)            #Asignar arrays a y,z
z = zeros(2,nt +1)
y[1] = 0.0              #Posición Inicial
y[2] = 2.0              #Angulo inicial
z[:,1] = y              #Ni idea
ni = 1                  #Contador

for t=0:dt:tfin         #Esto podria hacerse con un 't in range(n)' creo
    global ni = ni + 1 
    z[:,ni] = RK4(y,t,dt,q,b,w0)
    y[:] = z[:,ni]
end

plot(z[1,:],z[2,:],"r.")
axis([-pi,pi,-3,3])
xlabel(L"\theta")       #Aqui la 'L' entra en modo LaTeX
ylabel(L"\omega")
title("Espacio de fase")
savefig("Graf4.svg")
clf()

#Algo que no se que es, parece que hace arrays
#collect(1:8) .- 8/2 .-1
power = fft(z[2,1:nt]) |> fftshift;
freq = fftfreq(nt,2pi/(dt*w0)) |> fftshift;

plot(freq,log10.(abs2.(power/nt)))
xlabel(L"\omega/\omega_0")
ylabel(L"Power")
axis([0,6,-7,0])
savefig("Graf5.svg")

print("Todo Bien \n")