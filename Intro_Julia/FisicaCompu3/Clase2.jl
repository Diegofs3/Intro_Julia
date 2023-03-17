using PyPlot

function runge(x)
    return 1/(1 .+ 25*x.^2)
end

N = 100
x = range(-1,1,N)
xp = range(-1,1,10)

y = runge.(x)
yp = runge.(xp)
plot(x,y)
scatter(xp,yp, color = "red")
savefig("Graf1.png")
