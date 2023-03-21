######## Interpolacion de Lagrange ######
#Paquetes a usar#
using Pkg
#Pkg.add("Polynomials")
#Pkg.add("SpecialPolynomials")
#Pkg.add("DataInterpolations") Para la interpolacion en una variable
using Polynomials,SpecialPolynomials

#Primero trataremos de crear una función propia para la interpolación
#de Lagrange

function PolLagr(Datax,PNuevos)
    it = len(Datax) #Numero de iteraciones
    res = array(len(PNuevos))
    for i in range(it+1)
    #Pitatoria
        con = 1
        for j in range(it+1)
            if i ==j
                continue
            end
            con = con * (PNuevos -Datax[j])/(Datax[i]-Datax[j])
        end
        res[i] = con
    end
    return res
end
        


println("Todo bien")