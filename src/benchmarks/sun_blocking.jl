"""
    sun_blocking_pdgnep.jl

Sun-blocking differential game implementation using SpacecraftDynamics.jl
and DifferentialGames.jl.

# Problem Description
The inspector (pursuer) attempts to maintain line-of-sight (LoS) to a target
(evader) while the target maneuvers to block the LoS using the sun.

The reward is based on:
1. Evader-pursuer-sun angle (evader wants to maximize blocking)
2. Pursuer-evader distance (evader wants to stay at optimal viewing distance)

# Players
1. Evader (Target): Maneuvers to block sun while maintaining viewing distance
2. Pursuer (Inspector): Maneuvers to maintain LoS around evader

# Reference
Based on: https://github.com/mit-ll/spacegym-kspdg/tree/main/src/kspdg/sb1
"""

using SpacecraftDynamics
using DifferentialGames
using LinearAlgebra
using StaticArrays

# ============================================================================
# Problem Setup
# ============================================================================

"""
    create_sun_blocking_scenario(;
        altitude = 500e3,
        evader_mass = 150.0,
        pursuer_mass = 80.0,
        evader_thrust = 1.0,
        pursuer_thrust = 1.5,
        initial_separation = [500.0, 0.0, 0.0],
        pursuer_offset = [700.0, 100.0, 0.0]
    )

Create initial conditions for sun-blocking scenario.

# Keyword Arguments
- `altitude::Real`: Reference orbit altitude [m]
- `evader_mass::Real`: Evader spacecraft mass [kg]
- `pursuer_mass::Real`: Pursuer spacecraft mass [kg]
- `evader_thrust::Real`: Evader max thrust [N]
- `pursuer_thrust::Real`: Pursuer max thrust [N]
- `initial_separation::Vector`: Evader position relative to origin [m]
- `pursuer_offset::Vector`: Pursuer position relative to origin [m]

# Returns
Named tuple with scenario parameters
"""
function create_sun_blocking_scenario(;
    altitude::Real = 500e3,
    evader_mass::Real = 150.0,
    pursuer_mass::Real = 80.0,
    evader_thrust::Real = 1.0,
    pursuer_thrust::Real = 1.5,
    initial_separation::AbstractVector = [500.0, 0.0, 0.0],
    pursuer_offset::AbstractVector = [700.0, 100.0, 0.0]
)
    # Virtual chief orbit
    chief = VirtualChief(altitude=altitude)
    n = mean_motion(chief)
    
    # Create spacecraft (translational-only, nonlinear dynamics)
    # Evader
    evader_thrusters = ActuatorConfiguration(
        thrusters=default_research_thrusters(max_thrust=evader_thrust)
    )
    evader = Spacecraft(
        mass=evader_mass,
        inertia=Diagonal([10.0, 12.0, 8.0]),
        actuators=evader_thrusters,
        attitude_enabled=false,
        orbital_dynamics=:nonlinear_hcw
    )
    
    # Pursuer
    pursuer_thrusters = ActuatorConfiguration(
        thrusters=default_research_thrusters(max_thrust=pursuer_thrust)
    )
    pursuer = Spacecraft(
        mass=pursuer_mass,
        inertia=Diagonal([10.0, 12.0, 8.0]),
        actuators=pursuer_thrusters,
        attitude_enabled=false,
        orbital_dynamics=:nonlinear_hcw
    )
    
    # Initial states (6D: translational only)
    x0_evader = vcat(initial_separation, zeros(3))
    x0_pursuer = vcat(pursuer_offset, zeros(3))
    
    return (
        chief = chief,
        n = n,
        evader = evader,
        pursuer = pursuer,
        x0_evader = x0_evader,
        x0_pursuer = x0_pursuer,
        state_dims = [6, 6],
        control_dims = [6, 6]
    )
end

# ============================================================================
# Dynamics for PDGNEP
# ============================================================================

"""
    create_nonlinear_dynamics(spacecraft, n, μ)

Create dynamics function for nonlinear HCW with signature f(x, u, t).

# Arguments
- `spacecraft::Spacecraft`: Spacecraft configuration
- `n::Real`: Mean motion [rad/s]
- `μ::Real`: Gravitational parameter [m³/s²]

# Returns
Function with signature `(x, u, t) -> ẋ`
"""
function create_nonlinear_dynamics(spacecraft::Spacecraft, n::Real, μ::Real)
    # Create NonlinearHCWDynamics
    nonlin_dyn = NonlinearHCWDynamics(mass=spacecraft.mass, μ=μ)
    
    # Chief orbit parameters (circular)
    altitude = 500e3  # Default
    a = R_EARTH + altitude
    v_circ = n * a
    
    function dynamics(x::AbstractVector, u::AbstractVector, t::Real)
        # State: [r_rel (HCW), v_rel (HCW)] - 6D
        # Need to augment with chief state for nonlinear dynamics
        
        # Chief at current position (circular orbit in xy-plane)
        θ = n * t  # True anomaly for circular orbit
        r_chief = @SVector [a * cos(θ), a * sin(θ), 0.0]
        v_chief = @SVector [-v_circ * sin(θ), v_circ * cos(θ), 0.0]
        
        # Augmented state for nonlinear dynamics
        x_augmented = vcat(x, r_chief, v_chief)
        
        # Evaluate nonlinear dynamics
        ẋ_augmented = nonlin_dyn(x_augmented, u, t)
        
        # Extract only relative state derivative (first 6 elements)
        return ẋ_augmented[1:6]
    end
    
    return dynamics
end

# ============================================================================
# Objectives (Matching KSP-DG sb1_base.py::get_reward)
# ============================================================================

"""
    evader_stage_cost(x_evader, u_evader, x_pursuer, t;
                      target_viewing_distance, reward_decay_coef,
                      R, control_weight, sun_dir)

Evader's stage cost based on KSP-DG sb1 reward.

Maximizes (negative of):
    reward = -dot(u_evader_pursuer, u_sun_pursuer) * 
             exp(-reward_decay_coef * (distance - target_viewing_distance)^2)

Where:
- u_evader_pursuer: unit vector from pursuer to evader
- u_sun_pursuer: unit vector from pursuer to sun
- The dot product measures how well evader blocks the sun from pursuer's view

# Arguments
- `x_evader::AbstractVector`: Evader state [r; v]
- `u_evader::AbstractVector`: Evader control
- `x_pursuer::AbstractVector`: Pursuer state [r; v]
- `t::Real`: Time
- `target_viewing_distance::Real`: Optimal evader-pursuer distance [m]
- `reward_decay_coef::Real`: Distance penalty coefficient
- `R::AbstractMatrix`: Control cost matrix
- `control_weight::Real`: Control cost weight
- `sun_dir::AbstractVector`: Sun direction unit vector in HCW frame

# Returns
Cost value (evader minimizes negative reward)
"""
function evader_stage_cost(
    x_evader::AbstractVector,
    u_evader::AbstractVector,
    x_pursuer::AbstractVector,
    t::Real;
    target_viewing_distance::Real = 200.0,
    reward_decay_coef::Real = 1e-5,
    R::AbstractMatrix = Diagonal(0.01 * ones(6)),
    control_weight::Real = 0.1,
    sun_dir::AbstractVector = [1.0, 0.0, 0.0]  # Radial direction (sun)
)
    # Extract positions
    r_evader = x_evader[1:3]
    r_pursuer = x_pursuer[1:3]
    
    # Evader position relative to pursuer
    p_evader_pursuer = r_evader - r_pursuer
    d_evader_pursuer = norm(p_evader_pursuer)
    
    # Unit vectors
    if d_evader_pursuer > 1e-6
        u_evader_pursuer = p_evader_pursuer / d_evader_pursuer
    else
        u_evader_pursuer = [0.0, 0.0, 0.0]
    end
    
    u_sun_pursuer = sun_dir  # Sun direction (already unit vector)
    
    # Blocking reward (evader wants to maximize)
    # When dot product is -1: evader perfectly blocks sun
    # When dot product is +1: evader is opposite direction from sun
    blocking_reward = -dot(u_evader_pursuer, u_sun_pursuer)
    
    # Distance penalty (Gaussian around target viewing distance)
    distance_penalty = exp(-reward_decay_coef * (d_evader_pursuer - target_viewing_distance)^2)
    
    # Combined reward
    reward = blocking_reward * distance_penalty
    
    # Evader wants to MAXIMIZE reward, so MINIMIZE negative reward
    blocking_cost = -reward
    
    # Control cost
    control_cost = control_weight * dot(u_evader, R * u_evader)
    
    return blocking_cost + control_cost
end

"""
    pursuer_stage_cost(x_pursuer, u_pursuer, x_evader, t;
                       target_distance, Q_distance, R, control_weight)

Pursuer's stage cost: maintain desired distance while minimizing control.

Pursuer wants to:
1. Maintain optimal viewing distance
2. Minimize control effort

# Arguments
- `x_pursuer::AbstractVector`: Pursuer state [r; v]
- `u_pursuer::AbstractVector`: Pursuer control
- `x_evader::AbstractVector`: Evader state [r; v]
- `t::Real`: Time
- `target_distance::Real`: Desired pursuer-evader distance [m]
- `Q_distance::Real`: Distance tracking weight
- `R::AbstractMatrix`: Control cost matrix
- `control_weight::Real`: Control cost weight

# Returns
Cost value
"""
function pursuer_stage_cost(
    x_pursuer::AbstractVector,
    u_pursuer::AbstractVector,
    x_evader::AbstractVector,
    t::Real;
    target_distance::Real = 200.0,
    Q_distance::Real = 1.0,
    R::AbstractMatrix = Diagonal(0.01 * ones(6)),
    control_weight::Real = 0.1
)
    # Extract positions
    r_pursuer = x_pursuer[1:3]
    r_evader = x_evader[1:3]
    
    # Distance tracking
    distance = norm(r_evader - r_pursuer)
    distance_cost = Q_distance * (distance - target_distance)^2
    
    # Control cost
    control_cost = control_weight * dot(u_pursuer, R * u_pursuer)
    
    return distance_cost + control_cost
end

# ============================================================================
# Constraints
# ============================================================================

"""
    collision_avoidance_constraint(x, player_indices, params)

Shared constraint: maintain minimum separation between spacecraft.
"""
function collision_avoidance_constraint(
    x::AbstractVector,
    player_indices::Vector{Int},
    params::NamedTuple
)
    # Extract positions
    r_evader = x[1:3]
    r_pursuer = x[7:9]  # Offset by 6 (evader state dimension)
    
    distance = norm(r_pursuer - r_evader)
    min_distance = get(params, :min_distance, 10.0)
    
    # Constraint: distance >= min_distance
    # Reformulate as: min_distance² - distance² <= 0
    return min_distance^2 - distance^2
end

"""
    control_limits_constraint(u, params)

Private constraint: bound control inputs.
"""
function control_limits_constraint(
    u::AbstractVector,
    params::NamedTuple
)
    u_max = get(params, :u_max, 1.0)
    
    # Box constraints: -u_max <= u[i] <= u_max for each i
    # Reformulate as: |u[i]| - u_max <= 0
    return [abs(u[i]) - u_max for i in 1:length(u)]
end

# ============================================================================
# Main PDGNEP Construction
# ============================================================================

"""
    create_sun_blocking_pdgnep(;
        tf = 300.0,
        dt = 1.0,
        scenario = create_sun_blocking_scenario(),
        min_distance = 10.0,
        target_viewing_distance = 200.0,
        reward_decay_coef = 1e-5,
        R_evader = Diagonal(0.01 * ones(6)),
        R_pursuer = Diagonal(0.01 * ones(6)),
        control_weight = 0.1,
        Q_distance = 1.0
    )

Create sun-blocking PDGNEP matching KSP-DG sb1 formulation.

# Keyword Arguments
- `tf::Real`: Final time [s]
- `dt::Real`: Time step [s]
- `scenario`: Scenario from create_sun_blocking_scenario()
- `min_distance::Real`: Minimum separation constraint [m]
- `target_viewing_distance::Real`: Optimal viewing distance [m]
- `reward_decay_coef::Real`: Distance penalty coefficient (from KSP-DG)
- `R_evader::AbstractMatrix`: Evader control cost matrix
- `R_pursuer::AbstractMatrix`: Pursuer control cost matrix
- `control_weight::Real`: Overall control cost weight
- `Q_distance::Real`: Pursuer distance tracking weight

# Returns
`GameProblem` ready for solution with FALCON, CONDOR, etc.

# Notes
NO terminal costs - only stage costs (matching KSP-DG implementation)
"""
function create_sun_blocking_pdgnep(;
    tf::Real = 300.0,
    dt::Real = 1.0,
    scenario = create_sun_blocking_scenario(),
    min_distance::Real = 10.0,
    target_viewing_distance::Real = 200.0,
    reward_decay_coef::Real = 1e-5,
    R_evader::AbstractMatrix = Diagonal(0.01 * ones(6)),
    R_pursuer::AbstractMatrix = Diagonal(0.01 * ones(6)),
    control_weight::Real = 0.1,
    Q_distance::Real = 1.0
)
    # Extract scenario components
    evader = scenario.evader
    pursuer = scenario.pursuer
    x0_evader = scenario.x0_evader
    x0_pursuer = scenario.x0_pursuer
    n = scenario.n
    
    # Create dynamics functions
    dynamics_evader = create_nonlinear_dynamics(evader, n, μ_EARTH)
    dynamics_pursuer = create_nonlinear_dynamics(pursuer, n, μ_EARTH)
    
    # Create objectives with coupling (NO TERMINAL COSTS)
    # Evader's objective depends on pursuer state
    function evader_stage_cost_fn(x::AbstractVector, u::AbstractVector, t::Real)
        # Joint state: [x_evader; x_pursuer]
        x_evader = x[1:6]
        x_pursuer = x[7:12]
        u_evader = u[1:6]
        
        return evader_stage_cost(
            x_evader, u_evader, x_pursuer, t;
            target_viewing_distance=target_viewing_distance,
            reward_decay_coef=reward_decay_coef,
            R=R_evader,
            control_weight=control_weight
        )
    end
    
    # NO terminal cost for evader
    evader_terminal_cost_fn(x::AbstractVector) = 0.0
    
    # Pursuer's objective depends on evader state
    function pursuer_stage_cost_fn(x::AbstractVector, u::AbstractVector, t::Real)
        x_evader = x[1:6]
        x_pursuer = x[7:12]
        u_pursuer = u[7:12]
        
        return pursuer_stage_cost(
            x_pursuer, u_pursuer, x_evader, t;
            target_distance=target_viewing_distance,
            Q_distance=Q_distance,
            R=R_pursuer,
            control_weight=control_weight
        )
    end
    
    # NO terminal cost for pursuer
    pursuer_terminal_cost_fn(x::AbstractVector) = 0.0
    
    # Create Objective objects
    evader_objective = Objective(evader_stage_cost_fn, evader_terminal_cost_fn)
    pursuer_objective = Objective(pursuer_stage_cost_fn, pursuer_terminal_cost_fn)
    
    # Control limit constraints (private)
    evader_constraint_params = (u_max=1.0,)
    pursuer_constraint_params = (u_max=1.0,)
    
    evader_control_constraint = u -> control_limits_constraint(
        u, evader_constraint_params
    )
    pursuer_control_constraint = u -> control_limits_constraint(
        u, pursuer_constraint_params
    )
    
    # Create players
    player_evader = PlayerSpec(
        1,  # id
        6,  # state dimension
        6,  # control dimension
        x0_evader,
        dynamics_evader,
        evader_objective,
        [evader_control_constraint]  # private constraints
    )
    
    player_pursuer = PlayerSpec(
        2,  # id
        6,  # state dimension
        6,  # control dimension
        x0_pursuer,
        dynamics_pursuer,
        pursuer_objective,
        [pursuer_control_constraint]  # private constraints
    )
    
    players = [player_evader, player_pursuer]
    
    # Shared constraint: collision avoidance
    collision_params = (min_distance=min_distance,)
    collision_constraint_fn = (x, players) -> collision_avoidance_constraint(
        x, players, collision_params
    )
    
    shared_collision = SharedConstraint(
        collision_constraint_fn,
        [1, 2]  # Couples both players
    )
    
    # Create PDGNEP
    game = PDGNEProblem(
        players,
        [shared_collision],
        tf,
        dt
    )
    
    return game
end

# ============================================================================
# Example Usage
# ============================================================================

"""
Example: Create and solve sun-blocking game.
```julia
using SpacecraftDynamics
using DifferentialGames

# Create game problem
game = create_sun_blocking_pdgnep(
    tf=300.0,
    dt=1.0,
    min_distance=10.0,
    target_viewing_distance=200.0,
    reward_decay_coef=1e-5,
    control_weight=0.1
)

# Solve with FALCON (or other solver)
# solution = solve(game, FALCONSolver())

# Extract trajectories
# x_sol, u_sol = extract_solution(solution)
```
"""

# Export main functions
export create_sun_blocking_scenario
export create_sun_blocking_pdgnep
export evader_stage_cost, pursuer_stage_cost
export collision_avoidance_constraint