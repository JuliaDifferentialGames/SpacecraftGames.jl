# Lady-Guard-Bandit (LGB) three-player differential game.
#
# Recreates the KSPDG lbg1 scenario using realistic Earth-orbit dynamics
# (NonlinearHCW) instead of the KSP physics engine.
#
# Players
# -------
#   1 – Bandit  : Pursues the Lady while evading the Guard.
#   2 – Lady    : Evades the Bandit.
#   3 – Guard   : Pursues the Bandit to protect the Lady.
#
# Dynamics
# --------
#   Per-player separable NonlinearHCW in a 500 km LEO (Earth parameters).
#   State: 6-D [x, y, z, ẋ, ẏ, ż] in the rotating HCW/LVLH frame
#          x = radial-out, y = along-track (prograde), z = cross-track
#   Control: 3-D force [Fx, Fy, Fz] in HCW frame [N]
#
# Joint layout (18-D state, 9-D control)
# ----------------------------------------
#   state   : [x_bandit(1:6),  x_lady(7:12),  x_guard(13:18)]
#   control : [u_bandit(1:3),  u_lady(4:6),   u_guard(7:9)]
#   offsets (0-based): bandit=0, lady=6, guard=12
#
# Cost philosophy (per-step, all players minimise)
# -------------------------------------------------
#   Bandit : (d_BL / d_ref)²  +  w_bg / (d_BG / d_ref + 1)  +  R·‖u‖²/F²
#            ↑ approach Lady    ↑ barrier away from Guard
#   Lady   : exp(−d_BL² / σ²)  +  R·‖u‖²/F²
#            ↑ exponential penalty when Bandit is near
#   Guard  : (d_BG / d_ref)²   +  R·‖u‖²/F²
#            ↑ approach Bandit
#
# Initial conditions (Earth parameters, 500 km LEO)
# --------------------------------------------------
#   I1: Lady at origin, Guard 600 m prograde, Bandit on CW ellipse
#       (radial amplitude 500 m, along-track centre 1000 m prograde)
#       that passes through Lady's location after ~T/4 ≈ 23 min.
#   I2: All in the same circular orbit at different along-track offsets;
#       Bandit must thrust to close 2600 m while Guard tries to intercept.
#
# Reference: https://github.com/mit-ll/spacegym-kspdg/tree/main/src/kspdg/lbg1

# ===========================================================================
# Parameters
# ===========================================================================

"""
    LBGParameters

Physical and cost parameters for the Lady-Guard-Bandit game.

# Fields
- `altitude`        : Reference orbit altitude [m]
- `mass_bandit`     : Bandit mass [kg]
- `mass_lady`       : Lady mass [kg]
- `mass_guard`      : Guard mass [kg]
- `d_ref`           : Reference distance for cost normalisation [m]
- `w_bg`            : Bandit's Guard-avoidance barrier weight
- `sigma_lb`        : Lady-evasion Gaussian width [m]
- `w_evade`         : Lady evasion cost weight
- `w_intercept`     : Guard interception cost weight
- `R_bandit`        : Bandit control cost weight (normalised by max_thrust²)
- `R_lady`          : Lady control cost weight
- `R_guard`         : Guard control cost weight
- `max_thrust`      : Per-axis thrust bound [N]
- `collision_d_min` : Minimum separation for collision avoidance [m]
"""
Base.@kwdef struct LBGParameters
    altitude::Float64        = 500e3
    mass_bandit::Float64     = 100.0
    mass_lady::Float64       = 100.0
    mass_guard::Float64      = 100.0
    d_ref::Float64           = 100.0      # 100 m normalisation scale
    w_bg::Float64            = 1000.0     # Guard-avoidance barrier weight
    sigma_lb::Float64        = 500.0      # Lady evasion width [m]
    w_evade::Float64         = 1.0
    w_intercept::Float64     = 1.0
    R_bandit::Float64        = 0.01
    R_lady::Float64          = 0.01
    R_guard::Float64         = 0.01
    max_thrust::Float64      = 10.0       # N, per axis
    collision_d_min::Float64 = 50.0       # m
end

# ===========================================================================
# Initial-condition variants
# ===========================================================================

"""
    LBGInitVariant

Initial-condition variant for the LBG game.
- `LBG_I1`: Bandit on a CW ellipse that passes near the Lady (~23 min horizon).
- `LBG_I2`: All spacecraft in the same circular orbit; Bandit must thrust.
"""
@enum LBGInitVariant begin
    LBG_I1
    LBG_I2
end

"""
    lbg_initial_conditions(variant, params) -> NamedTuple

Compute HCW initial states for all three players.

Returns `(x0_bandit, x0_lady, x0_guard, x0_joint, n, a)`.
"""
function lbg_initial_conditions(variant::LBGInitVariant, params::LBGParameters)
    a = R_EARTH + params.altitude
    n = hcw_mean_motion(a; μ=μ_EARTH)

    if variant == LBG_I1
        # Lady at the HCW origin (chief's circular orbit).
        x0_lady  = zeros(6)
        # Guard 600 m prograde of Lady (y = along-track direction).
        x0_guard = [0.0, 600.0, 0.0, 0.0, 0.0, 0.0]
        # Bandit on a zero-drift CW ellipse:
        #   x(t) = A_r·cos(n·t),  y(t) = y_c − 2·A_r·sin(n·t)
        # At t=0: x=A_r, ẏ=-2·n·A_r.  Centre at y_c = 2·A_r so that the
        # ellipse passes through (0, 0) at t = T/4 ≈ 23 min.
        A_r = 500.0   # radial amplitude [m]
        y_c = 2.0 * A_r
        x0_bandit = [A_r, y_c, 0.0, 0.0, -2.0 * n * A_r, 0.0]

    elseif variant == LBG_I2
        # All three in the same circular orbit; Bandit must thrust to engage.
        x0_lady   = zeros(6)
        x0_guard  = [0.0, -600.0,  0.0, 0.0, 0.0, 0.0]   # 600 m retrograde
        x0_bandit = [0.0, -2600.0, 0.0, 0.0, 0.0, 0.0]   # 2600 m retrograde
    end

    return (
        x0_bandit = x0_bandit,
        x0_lady   = x0_lady,
        x0_guard  = x0_guard,
        x0_joint  = vcat(x0_bandit, x0_lady, x0_guard),
        n         = n,
        a         = a,
    )
end

# ===========================================================================
# Per-player NonlinearHCW dynamics (T4)
# ===========================================================================

# Returns f(xᵢ, uᵢ, p, t) → SVector{6}.  The chief orbit is computed
# analytically from the circular-orbit assumption at time t.
function _lbg_player_dynamics(mass::Float64, a::Float64, n::Float64)
    v_circ = n * a
    dyn    = NonlinearHCWDynamics(mass=mass, μ=μ_EARTH)

    function f(xᵢ, uᵢ, p, t)
        θ     = n * t
        r_c   = SVector{3,Float64}(a * cos(θ), a * sin(θ), 0.0)
        v_c   = SVector{3,Float64}(-v_circ * sin(θ), v_circ * cos(θ), 0.0)
        x_aug = vcat(SVector{6}(xᵢ), r_c, v_c)
        ẋ_aug = dyn(x_aug, SVector{3}(uᵢ), t)
        return SVector{6}(ẋ_aug[1], ẋ_aug[2], ẋ_aug[3],
                          ẋ_aug[4], ẋ_aug[5], ẋ_aug[6])
    end

    return f
end

# ===========================================================================
# Stage cost functions (T5)
# All receive the joint (x, u) and ignore the solver's p parameter.
# ===========================================================================

# Convenience: squared norm helper
@inline _sq(v) = dot(v, v)

function _lbg_bandit_stage_cost(x, u, params::LBGParameters, t)
    r_B = SVector{3}(x[1],  x[2],  x[3])
    r_L = SVector{3}(x[7],  x[8],  x[9])
    r_G = SVector{3}(x[13], x[14], x[15])
    u_B = SVector{3}(u[1],  u[2],  u[3])

    d_ref  = params.d_ref
    d_BL   = norm(r_B - r_L)
    d_BG   = norm(r_B - r_G)

    # Approach Lady: quadratic distance cost, normalised
    approach_cost  = (d_BL / d_ref)^2

    # Avoid Guard: soft barrier, grows as Guard closes in
    avoid_cost     = params.w_bg / (d_BG / d_ref + 1.0)

    # Control cost, normalised by max thrust²
    ctrl_cost      = params.R_bandit * _sq(u_B) / params.max_thrust^2

    return approach_cost + avoid_cost + ctrl_cost
end

function _lbg_lady_stage_cost(x, u, params::LBGParameters, t)
    r_B = SVector{3}(x[1], x[2], x[3])
    r_L = SVector{3}(x[7], x[8], x[9])
    u_L = SVector{3}(u[4], u[5], u[6])

    d_BL = norm(r_L - r_B)

    # Exponential penalty: large when Bandit is near, decays with distance
    evade_cost = params.w_evade * exp(-d_BL^2 / params.sigma_lb^2)
    ctrl_cost  = params.R_lady * _sq(u_L) / params.max_thrust^2

    return evade_cost + ctrl_cost
end

function _lbg_guard_stage_cost(x, u, params::LBGParameters, t)
    r_B = SVector{3}(x[1],  x[2],  x[3])
    r_G = SVector{3}(x[13], x[14], x[15])
    u_G = SVector{3}(u[7],  u[8],  u[9])

    d_ref  = params.d_ref
    d_BG   = norm(r_G - r_B)

    # Approach Bandit: quadratic distance cost, normalised
    intercept_cost = params.w_intercept * (d_BG / d_ref)^2
    ctrl_cost      = params.R_guard * _sq(u_G) / params.max_thrust^2

    return intercept_cost + ctrl_cost
end

# ===========================================================================
# Main game constructor (T7)
# ===========================================================================

"""
    create_lbg_game(; kwargs...) -> GameProblem{Float64}

Three-player Lady-Guard-Bandit differential game.

# Keyword Arguments
- `variant`  : `LBG_I1` or `LBG_I2` (default `LBG_I1`)
- `tf`, `dt` : Time horizon and step [s] (default 240 s, 1 s)
- `params`   : `LBGParameters` (default constructor)

# State / control layout
- 18-D joint state: [x_bandit(1:6), x_lady(7:12), x_guard(13:18)]
- 9-D joint control: [u_bandit(1:3), u_lady(4:6), u_guard(7:9)]

# Notes
The Bandit's cost depends on the Lady's and Guard's states simultaneously,
so all three objectives are non-separable (coupled) despite separable dynamics.
"""
function create_lbg_game(;
    variant::LBGInitVariant = LBG_I1,
    tf::Float64             = 240.0,
    dt::Float64             = 1.0,
    params::LBGParameters   = LBGParameters()
)
    ics = lbg_initial_conditions(variant, params)
    (; a, n, x0_bandit, x0_lady, x0_guard) = ics

    # --- Dynamics (captures params via closure) ---
    f_bandit = _lbg_player_dynamics(params.mass_bandit, a, n)
    f_lady   = _lbg_player_dynamics(params.mass_lady,   a, n)
    f_guard  = _lbg_player_dynamics(params.mass_guard,  a, n)

    # --- Objectives ---
    # Capture params in closures; ignore the solver-provided p argument.
    bandit_obj = PlayerObjective(1,
        NonlinearStageCost((x, u, p, t) -> _lbg_bandit_stage_cost(x, u, params, t)),
        NonlinearTerminalCost((x, p) -> 0.0)
    )
    lady_obj = PlayerObjective(2,
        NonlinearStageCost((x, u, p, t) -> _lbg_lady_stage_cost(x, u, params, t)),
        NonlinearTerminalCost((x, p) -> 0.0)
    )
    guard_obj = PlayerObjective(3,
        NonlinearStageCost((x, u, p, t) -> _lbg_guard_stage_cost(x, u, params, t)),
        NonlinearTerminalCost((x, p) -> 0.0)
    )

    # --- Control bounds (private per-player) ---
    F = params.max_thrust
    bandit_bounds = ControlBounds(1;
        control_offset = 0, control_dim = 3,
        lower = fill(-F, 3), upper = fill(F, 3))
    lady_bounds   = ControlBounds(2;
        control_offset = 3, control_dim = 3,
        lower = fill(-F, 3), upper = fill(F, 3))
    guard_bounds  = ControlBounds(3;
        control_offset = 6, control_dim = 3,
        lower = fill(-F, 3), upper = fill(F, 3))

    # --- Players ---
    player_bandit = Player{Float64}(1, 6, 3, x0_bandit, f_bandit, bandit_obj,
                                    [bandit_bounds])
    player_lady   = Player{Float64}(2, 6, 3, x0_lady,   f_lady,   lady_obj,
                                    [lady_bounds])
    player_guard  = Player{Float64}(3, 6, 3, x0_guard,  f_guard,  guard_obj,
                                    [guard_bounds])

    # --- Shared constraints: collision avoidance ---
    # State offsets (0-based): Bandit=0, Lady=6, Guard=12
    d_min = Float64(params.collision_d_min)
    col_BL = ProximityConstraint([1, 2];
        i_offset = 0, j_offset = 6,  pos_dim = 3, d_min = d_min)
    col_BG = ProximityConstraint([1, 3];
        i_offset = 0, j_offset = 12, pos_dim = 3, d_min = d_min)

    return DifferentialGame(
        [player_bandit, player_lady, player_guard], tf, dt;
        shared_constraints = [col_BL, col_BG]
    )
end

# ===========================================================================
# Bot policies (T8)
#
# These are fixed (non-game-theoretic) controllers for Lady and Guard used
# when simulating the game open-loop against a scripted opponent.
# Signature: policy(x_joint::Vector, t::Real, params::LBGParameters) → SVector{3}
# where the returned vector is the control force [Fx, Fy, Fz] in HCW [N].
# ===========================================================================

# ---------------------------------------------------------------------------
# LG0: Passive — Lady and Guard do not manoeuvre.
# ---------------------------------------------------------------------------

"""
    lbg_policy_lg0_lady(x_joint, t, params) -> SVector{3}

LG0 Lady policy: no thrust.
"""
lbg_policy_lg0_lady(x_joint, t, params::LBGParameters) = SVector(0.0, 0.0, 0.0)

"""
    lbg_policy_lg0_guard(x_joint, t, params) -> SVector{3}

LG0 Guard policy: no thrust.
"""
lbg_policy_lg0_guard(x_joint, t, params::LBGParameters) = SVector(0.0, 0.0, 0.0)

# ---------------------------------------------------------------------------
# LG1: Passive Lady, heuristic Guard pursuit (2-phase).
#
# Phase 1 (speed > threshold): zero out relative velocity (braking burn).
# Phase 2 (speed ≤ threshold): burn toward Bandit along line-of-sight.
#
# This is a stateless approximation of the KSPDG 3-phase policy; it omits
# the coasting phase because maintaining state across calls is outside scope.
# ---------------------------------------------------------------------------

"""
    lbg_policy_lg1_guard(x_joint, t, params) -> SVector{3}

LG1 Guard policy: 2-phase heuristic pursuit of the Bandit.

Phase 1: cancel relative velocity if closing speed is above threshold.
Phase 2: burn along the line-of-sight toward Bandit.
"""
function lbg_policy_lg1_guard(x_joint, t, params::LBGParameters)
    r_B = SVector{3}(x_joint[1],  x_joint[2],  x_joint[3])
    v_B = SVector{3}(x_joint[4],  x_joint[5],  x_joint[6])
    r_G = SVector{3}(x_joint[13], x_joint[14], x_joint[15])
    v_G = SVector{3}(x_joint[16], x_joint[17], x_joint[18])

    r_rel = r_B - r_G        # Bandit relative to Guard
    v_rel = v_B - v_G        # Relative velocity
    speed = norm(v_rel)

    min_speed_thresh = 0.2   # m/s (matches KSPDG default)
    F = params.max_thrust

    if speed > min_speed_thresh
        # Phase 1: decelerate relative motion
        u_dir = -v_rel / speed
        return SVector{3}(F .* u_dir)
    else
        d = norm(r_rel)
        if d > 1e-3
            # Phase 2: burn toward Bandit
            u_dir = r_rel / d
            return SVector{3}(F .* u_dir)
        else
            return SVector(0.0, 0.0, 0.0)
        end
    end
end

"""
    lbg_policy_lg1_lady(x_joint, t, params) -> SVector{3}

LG1 Lady policy: passive (same as LG0).
"""
lbg_policy_lg1_lady(x_joint, t, params::LBGParameters) = SVector(0.0, 0.0, 0.0)

# ---------------------------------------------------------------------------
# LG2: Evasive Lady (constant cross-track thrust), LG1 Guard.
#
# The Lady applies a constant out-of-plane (z-axis) burn.  In the HCW frame
# this generates oscillatory cross-track motion, making interception harder.
# ---------------------------------------------------------------------------

"""
    lbg_policy_lg2_lady(x_joint, t, params) -> SVector{3}

LG2 Lady policy: constant cross-track thrust to evade via out-of-plane motion.
"""
function lbg_policy_lg2_lady(x_joint, t, params::LBGParameters)
    return SVector(0.0, 0.0, params.max_thrust)
end

"""
    lbg_policy_lg2_guard(x_joint, t, params) -> SVector{3}

LG2 Guard policy: same 2-phase heuristic as LG1.
"""
lbg_policy_lg2_guard(x_joint, t, params::LBGParameters) =
    lbg_policy_lg1_guard(x_joint, t, params)
