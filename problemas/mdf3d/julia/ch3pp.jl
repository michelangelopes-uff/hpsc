# Compute Effective CONDUCTIVITY Properties
using JSON
using SparseArrays
using LinearAlgebra

# Model data struct:
struct Model
    nx::UInt64;
    ny::UInt64;
    nz::UInt64;
    voxelSize::Float64;
    refinement::UInt64;
    nMat::UInt16;
    rhsType::UInt8;
    solverType::UInt8;
    pcgTol::Float64;
    pcgIter::UInt64;
    matKeys::Vector{UInt16};
    matProp::Vector{Float64}; 
    nNodes::UInt64;
    nElems::UInt64;
    nDOFs::UInt64;
    DOFMap::Vector{UInt64}; 
    elemMatMap::Vector{UInt64}; 
    function Model(_nx, _ny, _nz, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap)
        new(_nx, _ny, _nz, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap);
    end
end

# Build Model:
function buildModel(_JsonFile::String, _RawFile::String)
    println("Construindo modelo")
    # Read Json file:
    nx, ny, nz, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp = readJSON(_JsonFile);
    # Read Raw file:
    elemMatMap = zeros(UInt64, nx * ny * nz * refinement * refinement * refinement);
    readRAW!(nx, ny, nz, refinement, matKeys, elemMatMap, _RawFile);
    # Update the parameters based on the given refinement level:
    nx *= refinement;
    ny *= refinement;
    nz *= refinement;
    nNodes = (nx + 1) * (ny + 1) * (nz + 1);
    nElems = (nx) * (ny) * (nz);
    DOFperNode = 1;
    nDOFs = nElems * DOFperNode;
    # Generate a map of Degree of Freedom:
    DOFMap = zeros(UInt64, nNodes);
    generateDOFMap!(nx, ny, nz, DOFMap);    
    # Build the Model:
    model = Model(nx, ny, nz, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp, nNodes, nElems, nDOFs, DOFMap, elemMatMap);
    return model
end

# Read JSON file:
function readJSON(_filename::String)
    # Open and read file:
    open(_filename, "r") do f
        data = JSON.parse(f);
        nx = data["image_dimensions"][1];
        ny = data["image_dimensions"][2];
        nz = data["image_dimensions"][3];
        refinement = 1;
        if haskey(data, "refinement");       refinement = data["refinement"];     end
        voxelSize = 1;
        if haskey(data, "voxel_size"); voxelSize = data["voxel_size"]; end
        rhsType = 0;
        if haskey(data, "type_of_rhs");      rhsType = data["type_of_rhs"];       end
        solverType = 0;
        if haskey(data, "type_of_solver");   solverType = data["type_of_solver"]; end
        pcgTol = 0.000001;
        if haskey(data, "solver_tolerance"); pcgTol = data["solver_tolerance"];   end       
        pcgIter = nx * ny * nz * refinement * refinement * refinement;
        if haskey(data, "number_of_iterations"); pcgIter = data["number_of_iterations"]; end
        nMat = data["number_of_materials"];
        materials = data["properties_of_materials"];
        matKeys = zeros(UInt16, 256);
        matProp = zeros(Float64, 256);
        for i = 1:nMat
            matKeys[convert(UInt8, materials[i][1]) + 1] = i;
            matProp[convert(UInt8, materials[i][1]) + 1] = convert(Float64, materials[i][2]);
        end
        materials = nothing;
        data = nothing;
        return nx, ny, nz, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp
    end    
end

# Read RAW file:
function readRAW!(_nx::Int, _ny::Int, _nz::Int, _refinement::Int, _matKeys::Vector{UInt16}, _elemMatMap::Vector{UInt64}, _filename::String)
    # Initializations
    rows = _ny; cols = _nx; slices = _nz;
    nelem = _nx * _ny; e = 0;
    line = _ny * _refinement * _refinement - _ny; buffer = 0;
    nx = _nx * _refinement;
    ny = _ny * _refinement;
    # Open and read file
    open(_filename, "r") do io
        bin_array = read(io);
        # Build the element material map based on the refinement level:
        for k = 1:slices
            for i = 1:rows
                for j = 1:cols
                    e = i + (j - 1) * _ny + (k - 1) * nelem;
                    buffer = _matKeys[bin_array[e] + 1];
                    for kk = (_refinement * (k - 1) + 1):(_refinement * k)
                        for ii = (_refinement * (i - 1) + 1):(_refinement * i)
                            for jj = (_refinement * (j - 1) + 1):(_refinement * j)
                                _elemMatMap[ii + (jj - 1) * ny + (kk - 1) * nx * ny] = buffer;
                            end
                        end
                    end          
                end
            end
        end
        bin_array = nothing;
    end
end

# Generate the Degree of Freedom Map:
function generateDOFMap!(_nx::Int, _ny::Int, _nz::Int, _DOFMap::Vector{UInt64})
    # Number the DOFs following the nodes from top to bottom and left to right:
    nElemS = _nx * _ny; 
    nNodeS = (_nx + 1) * (_ny + 1);
    @fastmath @inbounds @simd for n = 1:(nNodeS * (_nz + 1))
        i = (n - 1) % nNodeS;
        _DOFMap[n] = (i - (i ÷ (_ny + 1)) - _ny * ((i % (_ny + 1)) ÷ _ny)) % nElemS + (((n - 1) ÷ nNodeS) % _nz) * nElemS + 1;
    end
end

# Estimate memory consuption:
function estimateMemory(_model::Model)
    # elemMatMap = 16 bits * nElems
    # DOFMap = 64 bits * nNodes
    # RHS = 64 bits * nDOFs
    # PCGM Solver   / _solver == 0 / M d x q = 4 * 64 bits * _numDOFs / RHS = 64 bits * m_numDOFs
    # Direct Solver / _solver == 1 / K = 41 * 64 bits * _numElems (rough sparse estimative) / RHS = 3*64 bits * m_numDOFs
    mem = 0;
    if (_model.solverType == 0)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 5 * 64 * _model.nDOFs) / 8 / 1_000_000;
    elseif (_model.solverType == 1)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 3 * 64 * _model.nDOFs + 41 * 64 * _model.nElems) / 8 / 1_000_000;
    end
end

# Compute the element conductivity matrix for each material:
function elementConductivityMatrices!(_model::Model, _K::Array{Float64,3}, _B::Array{Float64,3})
    # Compute the matrices for each material:
    i = 0;
    for elemProps in _model.matProp[_model.matKeys .!= 0]
        i += 1;
        _K[:,:,i], _B[:,:,i] = C3D8ElementConductivity(elemProps);
    end
end

# Element C3D8 Conductivity - FEM
function C3D8ElementConductivity(_elemProps::Float64)
    # Initializations
    K = zeros(Float64, 8, 8);
    B = zeros(Float64, 3, 8);
    # Analytical K and B for voxel
    a = 1 / 3 * _elemProps; b = 1 / 12 * _elemProps;
    q = 1 / 4 * _elemProps;
    K[1,1] = +a; K[1,3] = -b; K[1,6] = -b; K[1,7] = -b; K[1,8] = -b;
    K[2,2] = +a; K[2,4] = -b; K[2,5] = -b; K[2,7] = -b; K[2,8] = -b;
    K[3,3] = +a; K[3,1] = -b; K[3,5] = -b; K[3,6] = -b; K[3,8] = -b;
    K[4,4] = +a; K[4,2] = -b; K[4,5] = -b; K[4,6] = -b; K[4,7] = -b;
    K[5,5] = +a; K[5,2] = -b; K[5,3] = -b; K[5,4] = -b; K[5,7] = -b;
    K[6,6] = +a; K[6,1] = -b; K[6,3] = -b; K[6,4] = -b; K[6,8] = -b;
    K[7,7] = +a; K[7,1] = -b; K[7,2] = -b; K[7,4] = -b; K[7,5] = -b;
    K[8,8] = +a; K[8,1] = -b; K[8,2] = -b; K[8,3] = -b; K[8,6] = -b;
    B[1,1] = -q; B[1,2] = +q; B[1,3] = +q; B[1,4] = -q;
    B[1,5] = -q; B[1,6] = +q; B[1,7] = +q; B[1,8] = -q;
    B[2,1] = -q; B[2,2] = -q; B[2,3] = +q; B[2,4] = +q;
    B[2,5] = -q; B[2,6] = -q; B[2,7] = +q; B[2,8] = +q;
    B[3,1] = -q; B[3,2] = -q; B[3,3] = -q; B[3,4] = -q;
    B[3,5] = +q; B[3,6] = +q; B[3,7] = +q; B[3,8] = +q;
    return K, B
end

# Compute the RHS: Boundary or Domain, _rhsType: Boundary = 0 || Domain = 1, _axis 0 = X || _axis 1 = Y || _axis 2 = Z
function computeRHS!(_model::Model, _RHS::Vector{Float64}, _axis::Int, _K::Array{Float64,3}, _B::Array{Float64,3})
    # Initializations
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; 
    N5 = 0; N6 = 0; N7 = 0; N8 = 0;
    nElemS = _model.nx * _model.ny; 
    nNodeS = (_model.nx + 1) * (_model.ny + 1);
    # Compute each RHS (_axis) based on boundary or domain data (_rhsType)
    if _model.rhsType == 1     # Boundary
        e = 0; 
        if _axis == 0  # _axis 0 = X
            deltaT = _model.nx;
            c = _model.nx;
            for dz = 1:_model.nz
                for r = 1:_model.ny
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
                    _RHS[_model.DOFMap[N1]] -= (_K[1,2,_model.elemMatMap[e]] + _K[1,3,_model.elemMatMap[e]] + _K[1,6,_model.elemMatMap[e]] + _K[1,7,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N2]] -= (_K[2,2,_model.elemMatMap[e]] + _K[2,3,_model.elemMatMap[e]] + _K[2,6,_model.elemMatMap[e]] + _K[2,7,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N3]] -= (_K[3,2,_model.elemMatMap[e]] + _K[3,3,_model.elemMatMap[e]] + _K[3,6,_model.elemMatMap[e]] + _K[3,7,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N4]] -= (_K[4,2,_model.elemMatMap[e]] + _K[4,3,_model.elemMatMap[e]] + _K[4,6,_model.elemMatMap[e]] + _K[4,7,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N5]] -= (_K[5,2,_model.elemMatMap[e]] + _K[5,3,_model.elemMatMap[e]] + _K[5,6,_model.elemMatMap[e]] + _K[5,7,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N6]] -= (_K[6,2,_model.elemMatMap[e]] + _K[6,3,_model.elemMatMap[e]] + _K[6,6,_model.elemMatMap[e]] + _K[6,7,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N7]] -= (_K[7,2,_model.elemMatMap[e]] + _K[7,3,_model.elemMatMap[e]] + _K[7,6,_model.elemMatMap[e]] + _K[7,7,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N8]] -= (_K[8,2,_model.elemMatMap[e]] + _K[8,3,_model.elemMatMap[e]] + _K[8,6,_model.elemMatMap[e]] + _K[8,7,_model.elemMatMap[e]]) * deltaT;
                end
            end
        elseif _axis == 1 # _axis 1 = Y
            deltaT = _model.ny; 
            r = 1; 
            for dz = 1:_model.nz
                for c = 1:_model.nx
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
                    _RHS[_model.DOFMap[N1]] -= (_K[1,3,_model.elemMatMap[e]] + _K[1,4,_model.elemMatMap[e]] + _K[1,7,_model.elemMatMap[e]] + _K[1,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N2]] -= (_K[2,3,_model.elemMatMap[e]] + _K[2,4,_model.elemMatMap[e]] + _K[2,7,_model.elemMatMap[e]] + _K[2,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N3]] -= (_K[3,3,_model.elemMatMap[e]] + _K[3,4,_model.elemMatMap[e]] + _K[3,7,_model.elemMatMap[e]] + _K[3,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N4]] -= (_K[4,3,_model.elemMatMap[e]] + _K[4,4,_model.elemMatMap[e]] + _K[4,7,_model.elemMatMap[e]] + _K[4,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N5]] -= (_K[5,3,_model.elemMatMap[e]] + _K[5,4,_model.elemMatMap[e]] + _K[5,7,_model.elemMatMap[e]] + _K[5,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N6]] -= (_K[6,3,_model.elemMatMap[e]] + _K[6,4,_model.elemMatMap[e]] + _K[6,7,_model.elemMatMap[e]] + _K[6,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N7]] -= (_K[7,3,_model.elemMatMap[e]] + _K[7,4,_model.elemMatMap[e]] + _K[7,7,_model.elemMatMap[e]] + _K[7,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N8]] -= (_K[8,3,_model.elemMatMap[e]] + _K[8,4,_model.elemMatMap[e]] + _K[8,7,_model.elemMatMap[e]] + _K[8,8,_model.elemMatMap[e]]) * deltaT;
                end
            end            
        elseif _axis == 2 # _axis 2 = Z
            deltaT = _model.nz;
            dz = _model.nz;
            for r = 1:_model.ny
                for c = 1:_model.nx
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;            
                    _RHS[_model.DOFMap[N1]] -= (_K[1,5,_model.elemMatMap[e]] + _K[1,6,_model.elemMatMap[e]] + _K[1,7,_model.elemMatMap[e]] + _K[1,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N2]] -= (_K[2,5,_model.elemMatMap[e]] + _K[2,6,_model.elemMatMap[e]] + _K[2,7,_model.elemMatMap[e]] + _K[2,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N3]] -= (_K[3,5,_model.elemMatMap[e]] + _K[3,6,_model.elemMatMap[e]] + _K[3,7,_model.elemMatMap[e]] + _K[3,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N4]] -= (_K[4,5,_model.elemMatMap[e]] + _K[4,6,_model.elemMatMap[e]] + _K[4,7,_model.elemMatMap[e]] + _K[4,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N5]] -= (_K[5,5,_model.elemMatMap[e]] + _K[5,6,_model.elemMatMap[e]] + _K[5,7,_model.elemMatMap[e]] + _K[5,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N6]] -= (_K[6,5,_model.elemMatMap[e]] + _K[6,6,_model.elemMatMap[e]] + _K[6,7,_model.elemMatMap[e]] + _K[6,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N7]] -= (_K[7,5,_model.elemMatMap[e]] + _K[7,6,_model.elemMatMap[e]] + _K[7,7,_model.elemMatMap[e]] + _K[7,8,_model.elemMatMap[e]]) * deltaT;
                    _RHS[_model.DOFMap[N8]] -= (_K[8,5,_model.elemMatMap[e]] + _K[8,6,_model.elemMatMap[e]] + _K[8,7,_model.elemMatMap[e]] + _K[8,8,_model.elemMatMap[e]]) * deltaT;
                end
            end
        end
    elseif _model.rhsType == 0  # Domain        
        for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;     
            _RHS[_model.DOFMap[N1]] += _B[_axis + 1,1,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N2]] += _B[_axis + 1,2,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N3]] += _B[_axis + 1,3,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N4]] += _B[_axis + 1,4,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N5]] += _B[_axis + 1,5,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N6]] += _B[_axis + 1,6,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N7]] += _B[_axis + 1,7,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N8]] += _B[_axis + 1,8,_model.elemMatMap[e]];
        end 
    end
end

# Direct Solver: [K] 64 bits * m_nDOFs * m_nDOFs 
function directMethod!(_model::Model, _x1::Vector{Float64}, _x2::Vector{Float64}, _x3::Vector{Float64}, _RHS1::Vector{Float64}, _RHS2::Vector{Float64}, _RHS3::Vector{Float64}, _K::Array{Float64,3})
    # Initializations
    K = spzeros(_model.nDOFs, _model.nDOFs);
    pElemDOFNum = zeros(UInt64, 8);
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; 
    N5 = 0; N6 = 0; N7 = 0; N8 = 0;
    nElemS = _model.nx * _model.ny; 
    nNodeS = (_model.nx + 1) * (_model.ny + 1);
    # Assembly system matrix:
    for e = 1:_model.nElems
        N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
        N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
        N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS; 
        pElemDOFNum[1] = _model.DOFMap[N1]; pElemDOFNum[2] = _model.DOFMap[N2]; 
        pElemDOFNum[3] = _model.DOFMap[N3]; pElemDOFNum[4] = _model.DOFMap[N4];
        pElemDOFNum[5] = _model.DOFMap[N5]; pElemDOFNum[6] = _model.DOFMap[N6]; 
        pElemDOFNum[7] = _model.DOFMap[N7]; pElemDOFNum[8] = _model.DOFMap[N8];        
        for i = 1:8
            for j = 1:8
                K[pElemDOFNum[i],pElemDOFNum[j]] += _K[i,j,_model.elemMatMap[e]];
            end
        end
    end
    # Solve for three rhs:
    _x1 .= K \ _RHS1;
    _x2 .= K \ _RHS2;
    _x3 .= K \ _RHS3;
end

# Jacobi Preconditioner: assembly || M
function jacobiPrecond!(_model::Model, _M::Vector{Float64}, _K::Array{Float64,3})
    # Initializations:
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    N5 = 0; N6 = 0; N7 = 0; N8 = 0;
    nElemS = _model.nx * _model.ny; 
    nNodeS = (_model.nx + 1) * (_model.ny + 1);
    for e = 1:_model.nElems
        N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
        N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
        N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;   
        _M[_model.DOFMap[N1]] += _K[1,1,_model.elemMatMap[e]];
        _M[_model.DOFMap[N2]] += _K[2,2,_model.elemMatMap[e]];
        _M[_model.DOFMap[N3]] += _K[3,3,_model.elemMatMap[e]];
        _M[_model.DOFMap[N4]] += _K[4,4,_model.elemMatMap[e]];
        _M[_model.DOFMap[N5]] += _K[5,5,_model.elemMatMap[e]];
        _M[_model.DOFMap[N6]] += _K[6,6,_model.elemMatMap[e]];
        _M[_model.DOFMap[N7]] += _K[7,7,_model.elemMatMap[e]];
        _M[_model.DOFMap[N8]] += _K[8,8,_model.elemMatMap[e]];
    end
    _M .= _M .\ 1;
end 

# Preconditioned Conjugate Gradient Method:
function pcg!(_model::Model, _x::Vector{Float64}, _r::Vector{Float64}, _M::Vector{Float64}, _K::Array{Float64,3})
    # Initializations
    d = zeros(Float64, _model.nDOFs);
    q = zeros(Float64, _model.nDOFs);
    pElemDOFNum = zeros(UInt64, 8);
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; 
    N5 = 0; N6 = 0; N7 = 0; N8 = 0;
    q_temp = 0;
    nElemS = _model.nx * _model.ny; 
    nNodeS = (_model.nx + 1) * (_model.ny + 1);
    # PCG Initialization:
    d .= _r;
    d .*= _M;
    delta_new = (_r' * d)[1,1];
    delta_0 = delta_new;
    i_max = _model.pcgIter;  
    ii = 0;
    # PCG Iterations:
    while (ii < i_max) && (abs(delta_new) > _model.pcgTol * _model.pcgTol * abs(delta_0)) # (maximum(abs.(_r))>_pcgTol)
        @fastmath @inbounds @simd for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
            pElemDOFNum[1] = _model.DOFMap[N1]; pElemDOFNum[2] = _model.DOFMap[N2]; 
            pElemDOFNum[3] = _model.DOFMap[N3]; pElemDOFNum[4] = _model.DOFMap[N4];
            pElemDOFNum[5] = _model.DOFMap[N5]; pElemDOFNum[6] = _model.DOFMap[N6]; 
            pElemDOFNum[7] = _model.DOFMap[N7]; pElemDOFNum[8] = _model.DOFMap[N8]; 
            for i = 1:8
                q_temp = 0;
                for j = 1:8
                    q_temp += _K[i,j,_model.elemMatMap[e]] * d[pElemDOFNum[j]];
                end
                q[pElemDOFNum[i]] += q_temp;
            end
        end
        alfa = delta_new / (d' * q)[1,1];
        d .*= alfa;
        _x .+= d;
        q .*= alfa;
        _r .-= q;
        q .= _r;
        q .*= _M;
        delta_old = delta_new;
        delta_new = (_r' * q)[1,1];
        beta = delta_new / delta_old;
        d .*= beta / alfa;
        d .+= q;
        q .*= 0;
        ii += 1;
    end
end

# Compute Flux-FEM Effective property
function femEffective(_model::Model, _T::Vector{Float64}, _axis::Int, _B::Array{Float64,3})
    # Initialization
    QX = 0; QY = 0; QZ = 0;
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; 
    N5 = 0; N6 = 0; N7 = 0; N8 = 0;
    nElemS = _model.nx * _model.ny; 
    nNodeS = (_model.nx + 1) * (_model.ny + 1);
    pElemDOFNum = zeros(UInt64, 8);
    bound = zeros(UInt8, 4);
    C = zeros(Float64, 3, 3);
    if _model.rhsType == 1  # Boundary
        if _axis == 0
            deltaT = _model.nx; aux = 0;
            bound[1] = 2; bound[2] = 3; bound[3] = 6; bound[4] = 7;
            for i = 1:_model.nz
                for eb = (nElemS + 1 - _model.ny + aux):(nElemS + aux)     
                    for j in bound
                        QX += (_B[1,j,_model.elemMatMap[eb]] * deltaT); 
                        QY += (_B[2,j,_model.elemMatMap[eb]] * deltaT); 
                        QZ += (_B[3,j,_model.elemMatMap[eb]] * deltaT); 
                    end
                end
                aux += nElemS;
            end 
        elseif _axis == 1
            deltaT = _model.ny;
            bound[1] = 3; bound[2] = 4; bound[3] = 7; bound[4] = 8;
            for eb = 1:(_model.ny):_model.nElems   
                for j in bound
                    QX += (_B[1,j,_model.elemMatMap[eb]] * deltaT); 
                    QY += (_B[2,j,_model.elemMatMap[eb]] * deltaT); 
                    QZ += (_B[3,j,_model.elemMatMap[eb]] * deltaT); 
                end
            end
        elseif _axis == 2
            deltaT = _model.nz;
            bound[1] = 5; bound[2] = 6; bound[3] = 7; bound[4] = 8;
            for eb = (_model.nElems - nElemS + 1):_model.nElems
                for j in bound
                    QX += (_B[1,j,_model.elemMatMap[eb]] * deltaT); 
                    QY += (_B[2,j,_model.elemMatMap[eb]] * deltaT); 
                    QZ += (_B[3,j,_model.elemMatMap[eb]] * deltaT); 
                end
            end  
        end  
        # Compute the effective properties for each test
        for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            pElemDOFNum[5] = N5; pElemDOFNum[6] = N6; pElemDOFNum[7] = N7; pElemDOFNum[8] = N8; 
            for i = 1:8
                QX += (_B[1,i,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[i]]]);
                QY += (_B[2,i,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[i]]]);
                QZ += (_B[3,i,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[i]]]);
            end
        end 
    elseif _model.rhsType == 0  # Domain
        t = zeros(Float64, 8);
        if (_axis == 0); t[2] = 1; t[3] = 1; t[6] = 1; t[7] = 1;
        elseif (_axis == 1); t[3] = 1; t[4] = 1; t[7] = 1; t[8] = 1;
        elseif (_axis == 2); t[5] = 1; t[6] = 1; t[7] = 1; t[8] = 1; end   
        for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            pElemDOFNum[5] = N5; pElemDOFNum[6] = N6; pElemDOFNum[7] = N7; pElemDOFNum[8] = N8; 
            for i = 1:8
                QX += (_B[1,i,_model.elemMatMap[e]] * (t[i] - _T[_model.DOFMap[pElemDOFNum[i]]]));
                QY += (_B[2,i,_model.elemMatMap[e]] * (t[i] - _T[_model.DOFMap[pElemDOFNum[i]]]));
                QZ += (_B[3,i,_model.elemMatMap[e]] * (t[i] - _T[_model.DOFMap[pElemDOFNum[i]]]));
            end
        end 
    end
    C[_axis + 1,1] = QX / _model.nElems; C[_axis + 1,2] = QY / _model.nElems; C[_axis + 1,3] = QZ / _model.nElems;
    return C
end

# -----------------
function mainCH(_arg)
    println("---------------------------", _arg);

    # Build the Model data struct:
    m_model = buildModel(_arg * ".json", _arg * ".raw");
    # Estimate Memory Consumption:
    estimateMemory(m_model);
    # SOLVE:
    # Compute the conductivity matrix for each Material:    
    m_K = zeros(Float64, 8, 8, m_model.nMat);
    m_B = zeros(Float64, 3, 8, m_model.nMat);
    elementConductivityMatrices!(m_model, m_K, m_B);
    if (m_model.solverType == 0) # Preconditioned Conjugate Gradient Method
        # Initialize the effective tensor, the right hand side, the inicial guess and the preconditioner:
        m_C = zeros(Float64, 3, 3);
        m_RHS = zeros(Float64, m_model.nDOFs);   
        m_X = zeros(Float64, m_model.nDOFs);
        m_M = zeros(Float64, m_model.nDOFs);
        # Compute the Jacobi preconditioner: 
        jacobiPrecond!(m_model, m_M, m_K);
        for axis = 0:2  # axis 0 = X || axis 1 = Y || axis 2 = Z
            # Compute the RHS: Boundary or Domain, m_rhsType: Boundary = 0 || Domain = 1, axis 0 = X || axis 1 = Y || axis 2 = Z
            computeRHS!(m_model, m_RHS, axis, m_K, m_B);            
            # Solver (to ensure optimal RAM usage we call GC before and after the PCGM):  
            GC.gc();            
            pcg!(m_model, m_X, m_RHS, m_M, m_K);
            GC.gc(); 
            # Compute Effective Property:
            m_C .+= femEffective(m_model, m_X, axis, m_B);
            m_RHS .*= 0;
            m_X .*= 0;
        end
    elseif (m_model.solverType == 1) # Direct Method
        # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y || axis 2 = Z
        m_RHS1 = zeros(Float64, m_model.nDOFs); m_RHS2 = zeros(Float64, m_model.nDOFs); m_RHS3 = zeros(Float64, m_model.nDOFs);
        computeRHS!(m_model, m_RHS1, 0, m_K, m_B);
        computeRHS!(m_model, m_RHS2, 1, m_K, m_B);
        computeRHS!(m_model, m_RHS2, 2, m_K, m_B);
        # Solver
        m_X1 = zeros(Float64, m_model.nDOFs); m_X2 = zeros(Float64, m_model.nDOFs); m_X3 = zeros(Float64, m_model.nDOFs);
        directMethod!(m_model, m_X1, m_X2, m_X3, m_RHS1, m_RHS2, m_RHS3, m_K);
        m_RHS1 = nothing; m_RHS2 = nothing; m_RHS3 = nothing;
        # Compute Effective Property:
        m_C = zeros(Float64, 3, 3);
        m_C .+= femEffective(m_model, m_X1, 0, m_B);
        m_C .+= femEffective(m_model, m_X2, 1, m_B);
        m_C .+= femEffective(m_model, m_X3, 2, m_B);
        m_X1 = nothing; m_X2 = nothing; m_X3 = nothing;
    end
    for i = 1:3
        # println("C[$i,:] = ", m_C[i,:]);
    end
    return m_model
end



# Starts application
# if length(ARGS) > 0
#     mainCH(ARGS[1])
# end
