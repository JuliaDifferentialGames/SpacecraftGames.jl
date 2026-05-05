# Formation control linear-quadratic game.
#
# Three spacecraft in a CW (Hill) reference frame cooperatively form a line
# along the x-axis.  Each agent has a 6-D HCW state [r; v] and a 3-D force
# control.  The shared linear dynamics are block-diagonal (separable per agent)
# but the costs couple all three agents, making this a genuine GNEP.
#
# Solvers: FNELQ (finite-horizon Nash LQ).

function _flq_blockdiag(blocks::AbstractVector{<:AbstractMatrix{T}}) where T
    nr = sum(size(b, 1) for b in blocks)
    nc = sum(size(b, 2) for b in blocks)
    out = zeros(T, nr, nc)
    r, c = 0, 0
    for b in blocks
        br, bc = size(b)
        out[r+1:r+br, c+1:c+bc] .= b
        r += br; c += bc
    end
    return out
end

"""
    create_formation_lq_game(; kwargs...) -> GameProblem{Float64}

Three-spacecraft CW formation-control game.

Each spacecraft controls its own 3-D force.  The shared state matrix is
block-diagonal HCW; each player's cost penalises its own position/velocity
error plus inter-agent coupling.

# Keyword Arguments
- `n_orbital`  : Mean motion [rad/s] (default 0.001)
- `tf`, `dt`   : Time horizon and step [s] (default 1000 s, 1 s)
- `Q_pos`      : Position tracking weight (default 1.0)
- `Q_vel`      : Velocity penalty weight (default 0.1)
- `Q_couple`   : Inter-agent coupling weight (default 0.5)
- `R_weight`   : Control cost weight (default 0.01)
- `r_des`      : Desired positions for each agent (default line along x-axis)
- `x0`         : Joint initial state (default triangle formation)
"""
function create_formation_lq_game(;
    n_orbital::Float64 = 0.001,
    tf::Float64        = 1000.0,
    dt::Float64        = 1.0,
    Q_pos::Float64     = 1.0,
    Q_vel::Float64     = 0.1,
    Q_couple::Float64  = 0.5,
    R_weight::Float64  = 0.01,
    r_des = [[-100.0, 0.0, 0.0], [0.0, 0.0, 0.0], [100.0, 0.0, 0.0]],
    x0 = vcat(
        [0.0, -50.0, 0.0, 0.0, 0.0, 0.0],
        [0.0,  50.0, 0.0, 0.0, 0.0, 0.0],
        [86.6,  0.0, 0.0, 0.0, 0.0, 0.0]
    )
)
    n_agents    = 3
    n_per_agent = 6
    n_total     = n_agents * n_per_agent

    # Shared HCW A matrix (block-diagonal)
    A_single = Matrix{Float64}(hcw_state_matrix(n_orbital))
    A = _flq_blockdiag([A_single for _ in 1:n_agents])

    # I₃ helper: 3×3 identity matrix as a concrete Matrix{Float64}.
    # (LinearAlgebra.I is UniformScaling — not callable as I(3).)
    eye3 = Matrix{Float64}(I, 3, 3)

    # Per-agent B matrices: force enters as acceleration in velocity rows
    B = [zeros(n_total, 3) for _ in 1:n_agents]
    for i in 1:n_agents
        off = (i - 1) * n_per_agent
        B[i][off+4:off+6, :] = eye3
    end

    # Per-agent Q matrices with inter-agent coupling.
    # Diagonal blocks use UniformScaling (+= s*I is valid for square views).
    # Off-diagonal blocks need a concrete matrix.
    Q = [zeros(n_total, n_total) for _ in 1:n_agents]
    for i in 1:n_agents
        oi = (i - 1) * n_per_agent
        Q[i][oi+1:oi+3, oi+1:oi+3] += Q_pos   * I     # UniformScaling ✓
        Q[i][oi+4:oi+6, oi+4:oi+6] += Q_vel   * I
        for j in 1:n_agents
            i == j && continue
            oj = (j - 1) * n_per_agent
            Q[i][oi+1:oi+3, oi+1:oi+3] +=  Q_couple * I
            Q[i][oj+1:oj+3, oj+1:oj+3] +=  Q_couple * I
            Q[i][oi+1:oi+3, oj+1:oj+3]  = -Q_couple * eye3   # concrete ✓
            Q[i][oj+1:oj+3, oi+1:oi+3]  = -Q_couple * eye3
        end
    end

    # Linear tracking terms (desired formation offsets)
    q = [zeros(n_total) for _ in 1:n_agents]
    for i in 1:n_agents
        oi = (i - 1) * n_per_agent
        q[i][oi+1:oi+3] = -2.0 * Q_pos * r_des[i]
    end

    R  = [diagm(fill(R_weight, 3)) for _ in 1:n_agents]
    Qf = [10.0 * Q[i] for i in 1:n_agents]

    return LQGameProblem(A, B, Q, R, Qf, Vector{Float64}(x0), tf; dt=dt, q=q)
end
