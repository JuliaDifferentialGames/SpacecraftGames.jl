module SpacecraftGames

using Reexport
using LinearAlgebra
using StaticArrays

# Re-export the full DifferentialGamesBase API so users only need
# `using SpacecraftGames` for everything.
@reexport using DifferentialGamesBase
@reexport using DifferentialGamesBaseSolvers

# SpacecraftDynamics is used internally by all game files.
# Not re-exported to avoid namespace pollution; import directly if needed.
using SpacecraftDynamics

# ---------------------------------------------------------------------------
# Visualization submodule
# ---------------------------------------------------------------------------
include("visualization/Visualization.jl")
using .Visualization
export Visualization

export SpacecraftGeometry, LightingConfig, CameraConfig, BackgroundConfig
export VideoExportSettings, GameVisualizationConfig
export cubesat_1u, cubesat_3u, cubesat_6u
export default_camera, default_lighting
export hd_video, uhd_video
export VizardScenario, write_vizard_file, visualize_with_vizard, run_vizard
export check_basilisk, install_basilisk, get_basilisk_info, remove_basilisk_venv

# ---------------------------------------------------------------------------
# Game benchmark: Formation control (3-agent LQ, CW dynamics)
# ---------------------------------------------------------------------------
include("benchmarks/formation_lq.jl")
export create_formation_lq_game

# ---------------------------------------------------------------------------
# Game benchmark: Sun-blocking pursuit-evasion (2-player, nonlinear HCW)
# ---------------------------------------------------------------------------
include("benchmarks/sun_blocking.jl")
export create_sun_blocking_scenario, create_sun_blocking_game

# ---------------------------------------------------------------------------
# Lady-Guard-Bandit game (3-player, nonlinear HCW)
# ---------------------------------------------------------------------------
include("games/lady_guard_bandit.jl")

# Parameters and variants
export LBGParameters, LBGInitVariant, LBG_I1, LBG_I2

# Game constructor and IC helper
export create_lbg_game, lbg_initial_conditions

# Bot policies (for open-loop simulation against scripted opponents)
export lbg_policy_lg0_lady,  lbg_policy_lg0_guard
export lbg_policy_lg1_lady,  lbg_policy_lg1_guard
export lbg_policy_lg2_lady,  lbg_policy_lg2_guard

end
