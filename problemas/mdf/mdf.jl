

function main(_file::String)
println("MDF")

conect = [
2   0   0   5;
3   1   0   6;
4   2   0   7;
0   3   0   8;
6   0   1   9;
7   5   2   10;
8   6   3   11;
0   7   4   12;
10  0   5   13;
11  9   6   14;
12  10  7   15;
0   11  8   16;
14  0   9   0;
15  13  10  0;
16  14  11  0;
0   15  12  0
]
cc = [
1   100;
1   75;
1   75;
1   0;
1   100;
0   0;
0   0;
1   0;
1   100;
0   0;
0   0;
1   0;
1   100;
1   25;
1   25;
1   0;
]
bloco = [4 -1 -1 -1 -1]

n,temp = size(conect)
A = zeros(Float64,n,n)
b = zeros(Float64,n,1)
# assembly
for i=1:n
    A[i,i] = bloco[1]
    for j=1:4
        col = conect[i,j]
        if col ≠ 0
            if cc[col,1] == 0
                A[i,col] = bloco[j+1]
            else
                b[i,1] -= bloco[j+1]*cc[col,2]
            end
        end
    end
end
# bc
for i=1:n
    if cc[i,1] == 1
        A[i,:] = zeros(Float64,1,n)
        A[i,i] = 1
        b[i,1] = cc[i,2]
    end
end

x = A\b
display(x)

end


if length(ARGS) == 1
    main(ARGS[1])
end