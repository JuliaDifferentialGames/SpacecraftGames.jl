module SpacecraftGames

# Core dependencies
using Reexport

# Base game formulations
@reexport using DifferentialGamesBase
@reexport using DifferentialGamesBaseSolvers

# Visualization submodule (includes everything)
include("visualization/Visualization.jl")

# Re-export visualization functionality
using .Visualization
export Visualization

# Re-export commonly used types and functions
export SpacecraftGeometry, LightingConfig, CameraConfig, BackgroundConfig
export VideoExportSettings, GameVisualizationConfig
export cubesat_1u, cubesat_3u, cubesat_6u
export default_camera, default_lighting
export hd_video, uhd_video
export VizardScenario, write_vizard_file, visualize_with_vizard, run_vizard
export check_basilisk, install_basilisk, get_basilisk_info, remove_basilisk_venv

end