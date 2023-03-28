function BSort(arr)
    n = length(arr)
    for i in 1:n-1
        for j in 1:n-1
            if arr[j] > arr[j+1]
                arr[j], arr[j+1] = arr[j+1], arr[j]
            end
        end
    end
    return arr
end

algo = rand(0:10,10)
print(algo)
println("\n")
res = BSort(algo)
print(res)
println("\n")