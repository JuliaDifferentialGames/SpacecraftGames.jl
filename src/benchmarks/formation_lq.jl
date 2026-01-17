using LinearAlgebra
using SparseArrays
using DifferentialGames
using DifferentialGamesBaseSolvers

function cwh_matrix(n_orbital::Float64)
    A = [
        0      0      0      1      0      0;
        0      0      0      0      1      0;
        0      0      0      0      0      1;
        3n_orbital^2  0      0      0      2n_orbital  0;
        0      0      0   -2n_orbital  0      0;
        0      0   -n_orbital^2  0      0      0
    ]
    return A
end

# Manual block diagonal construction for dense matrices
function manual_blockdiag(blocks::Vector{Matrix{T}}) where T
    n = sum(size(b, 1) for b in blocks)
    m = sum(size(b, 2) for b in blocks)
    result = zeros(T, n, m)
    
    row_offset = 0
    col_offset = 0
    for block in blocks
        nr, nc = size(block)
        result[row_offset+1:row_offset+nr, col_offset+1:col_offset+nc] = block
        row_offset += nr
        col_offset += nc
    end
    
    return result
end

# Parameters
n_orbital = 0.001  # rad/s
n_agents = 3
n_per_agent = 6
n_total = 18

# Dynamics
A_single = cwh_matrix(n_orbital)
A = manual_blockdiag([A_single for _ in 1:n_agents])

B = [zeros(n_total, 3) for _ in 1:n_agents]
for i in 1:n_agents
    offset = (i-1) * n_per_agent
    B[i][offset+4:offset+6, :] = Matrix{Float64}(I, 3, 3)
end

# Desired formation (line along x-axis)
r_des = [[-100.0, 0.0, 0.0], [0.0, 0.0, 0.0], [100.0, 0.0, 0.0]]

# Initial state (triangle)
x0 = vcat(
    [0.0, -50.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 50.0, 0.0, 0.0, 0.0, 0.0],
    [86.6, 0.0, 0.0, 0.0, 0.0, 0.0]
)

# Weights
Q_pos, Q_vel, Q_couple, R_weight = 1.0, 0.1, 0.5, 0.01

# Build Q matrices with coupling
Q = Vector{Matrix{Float64}}(undef, n_agents)
for i in 1:n_agents
    Q[i] = zeros(n_total, n_total)
    
    # Self-cost
    offset_i = (i-1) * n_per_agent
    Q[i][offset_i+1:offset_i+3, offset_i+1:offset_i+3] = Q_pos * I(3)
    Q[i][offset_i+4:offset_i+6, offset_i+4:offset_i+6] = Q_vel * I(3)
    
    # Coupling with all other agents
    for j in 1:n_agents
        if i != j
            offset_j = (j-1) * n_per_agent
            
            # ||ri - rj||² coupling
            Q[i][offset_i+1:offset_i+3, offset_i+1:offset_i+3] += Q_couple * I(3)
            Q[i][offset_j+1:offset_j+3, offset_j+1:offset_j+3] += Q_couple * I(3)
            Q[i][offset_i+1:offset_i+3, offset_j+1:offset_j+3] = -Q_couple * I(3)
            Q[i][offset_j+1:offset_j+3, offset_i+1:offset_i+3] = -Q_couple * I(3)
        end
    end
end

# Linear tracking terms
q = [zeros(n_total) for _ in 1:n_agents]
for i in 1:n_agents
    offset_i = (i-1) * n_per_agent
    q[i][offset_i+1:offset_i+3] = -2 * Q_pos * r_des[i]
end

# Control and terminal costs
R = [diagm(fill(R_weight, 3)) for _ in 1:n_agents]
Qf = [10.0 * Q[i] for i in 1:n_agents]

# Create game
game = LQGameProblem(A, B, Q, R, Qf, x0, 1000.0; dt=1.0, q=q)

# Solve via FNELQ
solver = FNELQ()
sol = solve(game, solver; verbose=true)
