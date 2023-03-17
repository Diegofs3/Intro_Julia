####### Paquetes a usar ######3
using Pkg
#Pkg.add("PyPlot")
#Pkg.add("Polynomials")
#Pkg.add("CubicSplines")
using PyPlot,Polynomials,CubicSplines


####Clase 2#####

function runge(x)
    return 1/(1 .+ 25*x.^2)
end

N = 40
x = range(-1,1,N)
xp = range(-1,1,40)

y = runge.(x)
yp = runge.(xp)
plot(x,y)
scatter(xp,yp, color = "red")
savefig("Graf1.svg")
clf() 
#### Clase 3 ######

f1 = fit(xp,yp,7)
yf1 = f1.(xp)
plot(xp,yf1)
scatter(xp,yp,color = "red")
savefig("Graf2.svg")
clf()

##### Spline ####

ajus = CubicSpline(x,y)

xs = range(x[1], stop =x[end], length = 1000)
ys = ajus[xs]
plot(x,y, "o")
plot(xs,ys)
savefig("Graf3.svg")
clf()
