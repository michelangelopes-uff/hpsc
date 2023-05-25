include("./ch3pp.jl")
using Printf
using SparseArrays
using LinearAlgebra
using JSON

ITERATION_LIMIT = 1000

function main(conect, cc, nx, ny, nz)
    # conect = [
        # 0   0   4   2   0   10  # 1
        # 1   0   5   3   0   11  # 2
        # 2   0   6   0   0   12  # 3
        # 0   1   7   5   0   13  # 4
        # 4   2   8   6   0   14  # 5
        # 5   3   9   0   0   15  # 6
        # 0   4   0   8   0   16  # 7
        # 7   5   0   9   0   17  # 8
        # 8   6   0   0   0   18  # 9
        # 0   0   13  11  1   19  # 10
        # 10  0   14  11  2   20  # 11
        # 11  0   15  0   3   21  # 12
        # 0   10  16  14  4   22  # 13
        # 13  11  17  15  5   23  # 14
        # 14  12  18  0   6   24  # 15
        # 0   13  0   17  7   25  # 16
        # 16  14  0   18  8   26  # 17
        # 17  15  0   0   9   27  # 18
        # 0   0   22  20  10  0   # 19
        # 19  0   23  21  11  0   # 20
        # 20  0   24  0   12  0   # 21
        # 0   19  25  23  13  0   # 22
        # 22  20  26  24  14  0   # 23
        # 23  21  27  0   15  0   # 24
        # 0   22  0   26  16  0   # 25
        # 25  23  0   27  17  0   # 26
        # 26  24  0   0   18  0   # 27
        # ]

    # cc = [
        # 1  100 # 1
        # 1  100 # 2
        # 1  100 # 3
        # 1  0   # 4
        # 1  0   # 5
        # 1  0   # 6
        # 1  0   # 7
        # 1  0   # 8
        # 1  0   # 9
        # 1  100 # 10
        # 1  100 # 11
        # 1  100 # 12
        # 1  0   # 13
        # 0  0   # 14
        # 1  0   # 15
        # 1  0   # 16
        # 1  0   # 17
        # 1  0   # 18
        # 1  100 # 19
        # 1  100 # 20
        # 1  100 # 21
        # 1  0   # 22
        # 1  0   # 23
        # 1  0   # 24
        # 1  0   # 25
        # 1  0   # 26
        # 1  0   # 27
    # ]

    nn = size(conect)[1]
    A = spzeros( nn, nn)
    # A = zeros(Float64, nn, nn) #usar sparce para guardar a matriz A
    b = zeros(Float64, nn, 1)
    bloco = [1, 1, 1, 1, 1, 1]
    println("criando matriz A")
    for e = 1:nn
        if(cc[e,1] == 0)
            A[e,e] = -6
            A[e,conect[e,1]] = 1
            A[e,conect[e,2]] = 1
            A[e,conect[e,3]] = 1
            A[e,conect[e,4]] = 1
            A[e,conect[e,5]] = 1
            A[e,conect[e,6]] = 1
        else
            A[e,e] = 1
            b[e,1] = cc[e,2]
        end
    end
    # @show(b)
    # x = A\b
    # A = sparse(A)
    # @show A
    println("calculando gaussseidel")
    x = gaussseidel(A, b)
    println("fim do calculo")
    # @show nx, ny , nz
    imagDim = nx * ny
    layers = (nz - 1) / 2
    posInicial = Int64((imagDim * layers)+1)
    posFinal = Int64(posInicial + imagDim - 1)
    imagem = zeros(imagDim)

    imagem = x[posInicial: posFinal]
    imagem = reshape(imagem, (nx,ny))
    # imagem = transpose(imagem)
    # @show val1 = convert(Int64,imagem[6] )
    # imagem = round(imagem)
    outputRes(imagem)
    # @show imagem
    # 
    # k = x[25,37]
    # imagem = reshape(k, (nx,ny))
    # @show k

end

function teste2(_model::Model, deltaT)
    cor = zeros(Int16, _model.nMat);
    count = 1
    for m = 1:256
        if (_model.matKeys[m] > 0)
            cor[count] = m-1
            count += 1
        end
    end
    if (cor[_model.nMat] == 0 )
        cor = circshift(cor,1)
    end
    nz = (convert(Int64,_model.nz)-1)/2
    @show nz
    dT = parse(Int64, deltaT)/(nz)
    return criaConnect(_model, cor, dT, nz)
end

# <: (e - 1) % x == 0 ---> invalido
# V: (e - 1) % (x * y) >= (x * y) - x ---> inválido
# ^: (e - 1) % (x * y) < x ---> inválido
# >: (e - 1) % x == (x - 1) ---> invalido
# s: r < 1 ---> invalido
# w: r > x*y*z ---> invalido

# 19 <: 19 - 1 = 18
# 19 V: 19 + x = 19 + 4 = 23
# 16 ^: 19 - x = 19 - 4 = 15
# 19 >: 19 + 1 = 20
# 19 s: 19 - x*y = 19 - 4*3 = 7
# 19 w: 19 + x*y = 19 + 4*3 = 31


function criaConnect(_model::Model, cor, deltaT, layers)
    layerAtual = layers
    @show layerAtual
    dim = convert(Int64,_model.nx * _model.ny * _model.nz)
    # @show dim
    xy = convert(Int64,_model.nx * _model.ny)
    x = convert(Int64,_model.nx)
    isContorno = false
    conncet = zeros(Int64, dim, 6)
    cc = zeros(Int64, dim, 2)
    indConnect = 1
    for e = 1:dim
        isContorno = false
        if (mod(e - 1, x) == 0) # direita
            isContorno = true
        else
            conncet[indConnect,1] = e -1
        end

        if (mod(e - 1, xy) >= xy - x) # baixo
            isContorno = true
        else
            conncet[indConnect,2] = e + x
        end

        if (mod((e - 1), xy) < x) # cima
            isContorno = true
        else
            conncet[indConnect,3] = e - x
        end

        if (mod((e - 1), x) == (x - 1)) # esquerda
            isContorno = true
        else
            conncet[indConnect,4] = e + 1
        end

        s = e - xy
        if (s < 1)
            isContorno = true
        else
            conncet[indConnect,5] = s
        end

        w = e + xy
        if (w > dim)
            isContorno = true
        else
            conncet[indConnect,6] = w
        end
        if (isContorno)
            if (mod(e-1,xy) == 0 && e != 1)
                layerAtual -= 1
            end
            cc[e,1] = 1
            cc[e,2] = cor[_model.elemMatMap[indConnect]] + trunc(deltaT) * layerAtual
        end
        indConnect += 1
    end
    return conncet, cc

end

function gaussseidel(A, b)
    # println("https://en.wikipedia.org/wiki/Gauss–Seidel_method")
    x = zeros(size(b)[1])
    for it_count = 1:ITERATION_LIMIT
        x_new = zeros(size(x)[1])
        for i = 1:size(A)[1]
            s1 = dot(A[i,1:i-1], x_new[1:i-1])
            s2 = dot(A[i, i+1:size(A)[2]], x[i+1:size(x)[1]])
            x_new[i] = (b[i] - s1 - s2) / A[i, i]
        end
        if (isapprox(x, x_new))
            break
        end
        x = x_new
    end
    return x
end

function outputRes(_res)
    dict = Dict()
    push!(dict,"imagem"=>_res)
    open("imagem.json","w") do f
        JSON.print(f,dict)
    end
end

if length(ARGS) > 1
    @time model = mainCH(ARGS[1])
    println("modelo construido")
    println("construindo matriz conncet")
    @time conncet, cc = teste2(model, ARGS[2])
    println("matriz construida")
    println("começando calculo")
    @time main(conncet, cc, convert(Int64,model.nx), convert(Int64,model.ny), convert(Int64,model.nz))

else
    println("passe o nome do arquivo sem extenção e o deltaT")

end