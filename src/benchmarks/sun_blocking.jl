# Sun-blocking pursuit-evasion differential game.
#
# Two spacecraft orbit Earth at 500 km LEO.  The Evader (Sun-blocker) tries to
# position itself between the Pursuer (Inspector) and the Sun.  The Pursuer
# tries to maintain a desired viewing distance from the Evader.
#
# Based on: spacegym-kspdg sb1_base.py
# Dynamics : NonlinearHCW, 6-D state [r_rel; v_rel] per player, 3-D force control.
# Players  : 1 = Evader, 2 = Pursuer

# ---------------------------------------------------------------------------
# Scenario builder
# ---------------------------------------------------------------------------

"""
    create_sun_blocking_scenario(; kwargs...) -> NamedTuple

Return orbit parameters and initial conditions for the sun-blocking game.

# Keyword Arguments
- `altitude`            : Reference orbit altitude [m] (default 500 km)
- `evader_mass`         : Evader mass [kg] (default 150)
- `pursuer_mass`        : Pursuer mass [kg] (default 80)
- `initial_separation`  : Evader position in HCW [m] (default [500, 0, 0])
- `pursuer_offset`      : Pursuer position in HCW [m] (default [700, 100, 0])
"""
function create_sun_blocking_scenario(;
    altitude::Float64           = 500e3,
    evader_mass::Float64        = 150.0,
    pursuer_mass::Float64       = 80.0,
    initial_separation::Vector{Float64} = [500.0, 0.0, 0.0],
    pursuer_offset::Vector{Float64}     = [700.0, 100.0, 0.0]
)
    a = R_EARTH + altitude
    n = hcw_mean_motion(a; μ=μ_EARTH)

    x0_evader  = vcat(initial_separation, zeros(3))
    x0_pursuer = vcat(pursuer_offset,     zeros(3))

    return (
        altitude     = altitude,
        a            = a,
        n            = n,
        evader_mass  = evader_mass,
        pursuer_mass = pursuer_mass,
        x0_evader    = x0_evader,
        x0_pursuer   = x0_pursuer,
    )
end

# ---------------------------------------------------------------------------
# Nonlinear HCW dynamics wrapper (3-D force input, 6-D relative state output)
# ---------------------------------------------------------------------------

function _sb_player_dynamics(mass::Float64, a::Float64, n::Float64)
    v_circ = n * a
    dyn    = NonlinearHCWDynamics(mass=mass, μ=μ_EARTH)

    # Signature required by PDGNEProblem: fᵢ(xᵢ, uᵢ, p, t) → ẋᵢ
    function f(xᵢ, uᵢ, p, t)
        θ     = n * t
        r_c   = SVector{3,Float64}(a * cos(θ), a * sin(θ), 0.0)
        v_c   = SVector{3,Float64}(-v_circ * sin(θ), v_circ * cos(θ), 0.0)
        x_aug = vcat(SVector{6}(xᵢ), r_c, v_c)       # 12-D augmented state
        ẋ_aug = dyn(x_aug, SVector{3}(uᵢ), t)
        return SVector{6}(ẋ_aug[1], ẋ_aug[2], ẋ_aug[3],
                          ẋ_aug[4], ẋ_aug[5], ẋ_aug[6])
    end

    return f
end

# ---------------------------------------------------------------------------
# Stage costs
#
# Joint state layout : [x_evader(6D), x_pursuer(6D)]  (indices 1:6, 7:12)
# Joint control layout: [u_evader(3D), u_pursuer(3D)]  (indices 1:3, 4:6)
# Sun direction (HCW) : radial-outward = x-axis = [1,0,0]
# ---------------------------------------------------------------------------

function _sb_evader_stage_cost(x, u, p, t;
    target_distance  = 200.0,
    reward_decay     = 1e-5,
    control_weight   = 0.1,
    sun_dir          = SVector(1.0, 0.0, 0.0)
)
    r_evader  = SVector{3}(x[1], x[2], x[3])
    r_pursuer = SVector{3}(x[7], x[8], x[9])
    u_evader  = SVector{3}(u[1], u[2], u[3])

    Δr = r_evader - r_pursuer
    d  = norm(Δr)

    # Smooth regularisation: avoids a type-unstable ternary during ForwardDiff
    # (the two branches would return SVector{3,Dual} vs SVector{3,Float64}).
    # 1e-9 m regularisation is negligible for any physically realistic separation.
    û_ep = Δr / (d + 1e-9)

    # Negative dot product → -1 when evader is exactly between Sun and Pursuer
    blocking_reward = -dot(û_ep, sun_dir)
    distance_factor = exp(-reward_decay * (d - target_distance)^2)

    # Evader minimises −reward
    blocking_cost = -(blocking_reward * distance_factor)
    ctrl_cost     = control_weight * dot(u_evader, u_evader)

    return blocking_cost + ctrl_cost
end

function _sb_pursuer_stage_cost(x, u, p, t;
    target_distance = 200.0,
    Q_distance      = 1.0,
    control_weight  = 0.1
)
    r_evader  = SVector{3}(x[1], x[2], x[3])
    r_pursuer = SVector{3}(x[7], x[8], x[9])
    u_pursuer = SVector{3}(u[4], u[5], u[6])

    d            = norm(r_evader - r_pursuer)
    distance_cost = Q_distance * (d - target_distance)^2
    ctrl_cost     = control_weight * dot(u_pursuer, u_pursuer)

    return distance_cost + ctrl_cost
end

# ---------------------------------------------------------------------------
# Main constructor
# ---------------------------------------------------------------------------

"""
    create_sun_blocking_game(; kwargs...) -> GameProblem{Float64}

Two-player sun-blocking differential game (nonlinear HCW, 3-D force control).

Player 1 = Evader: manoeuvres to block the Sun from Pursuer's line of sight.
Player 2 = Pursuer: maintains desired viewing distance from Evader.

# Keyword Arguments
- `tf`, `dt`             : Time horizon / step [s] (default 300 s, 1 s)
- `scenario`             : Scenario NamedTuple from `create_sun_blocking_scenario()`
- `min_distance`         : Collision-avoidance separation [m] (default 10)
- `target_distance`      : Desired Pursuer-Evader separation [m] (default 200)
- `reward_decay`         : Gaussian width for evader blocking reward [1/m²] (default 1e-5)
- `control_weight`       : Control cost weight for both players (default 0.1)
- `Q_distance`           : Pursuer distance-tracking weight (default 1.0)
- `max_thrust`           : Per-axis thrust bound [N] (default 1.0)
"""
function create_sun_blocking_game(;
    tf::Float64             = 300.0,
    dt::Float64             = 1.0,
    scenario                = create_sun_blocking_scenario(),
    min_distance::Float64   = 10.0,
    target_distance::Float64 = 200.0,
    reward_decay::Float64   = 1e-5,
    control_weight::Float64 = 0.1,
    Q_distance::Float64     = 1.0,
    max_thrust::Float64     = 1.0
)
    (; altitude, a, n, evader_mass, pursuer_mass, x0_evader, x0_pursuer) = scenario

    # --- Dynamics ---
    f_evader  = _sb_player_dynamics(evader_mass,  a, n)
    f_pursuer = _sb_player_dynamics(pursuer_mass, a, n)

    # --- Objectives ---
    # Capture cost params in closures; ignore p from the solver (unused here)
    evader_stage = NonlinearStageCost(
        (x, u, p, t) -> _sb_evader_stage_cost(x, u, p, t;
            target_distance = target_distance,
            reward_decay    = reward_decay,
            control_weight  = control_weight)
    )
    evader_term = NonlinearTerminalCost((x, p) -> 0.0)
    evader_obj  = PlayerObjective(1, evader_stage, evader_term)

    pursuer_stage = NonlinearStageCost(
        (x, u, p, t) -> _sb_pursuer_stage_cost(x, u, p, t;
            target_distance = target_distance,
            Q_distance      = Q_distance,
            control_weight  = control_weight)
    )
    pursuer_term = NonlinearTerminalCost((x, p) -> 0.0)
    pursuer_obj  = PlayerObjective(2, pursuer_stage, pursuer_term)

    # --- Control bounds (private constraints) ---
    evader_bounds  = ControlBounds(1;
        control_offset = 0, control_dim = 3,
        lower = fill(-max_thrust, 3), upper = fill(max_thrust, 3))
    pursuer_bounds = ControlBounds(2;
        control_offset = 3, control_dim = 3,
        lower = fill(-max_thrust, 3), upper = fill(max_thrust, 3))

    # --- Players ---
    player_evader  = Player{Float64}(1, 6, 3, x0_evader,  f_evader,  evader_obj,
                                     [evader_bounds])
    player_pursuer = Player{Float64}(2, 6, 3, x0_pursuer, f_pursuer, pursuer_obj,
                                     [pursuer_bounds])

    # --- Shared constraint: collision avoidance ---
    # State offsets (0-based): evader=0, pursuer=6
    collision = ProximityConstraint([1, 2];
        i_offset = 0, j_offset = 6, pos_dim = 3,
        d_min    = Float64(min_distance))

    return DifferentialGame(
        [player_evader, player_pursuer], tf, dt;
        shared_constraints = [collision]
    )
end
