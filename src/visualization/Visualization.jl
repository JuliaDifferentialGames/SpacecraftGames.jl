# src/visualization/Visualization.jl

"""
Main visualization module for SpacecraftGames.jl
"""
module Visualization

# Re-export core types and functions
include("VisualizationCore.jl")
using .VisualizationCore

# Basilisk setup (include first, before VizardInterface)
include("BasiliskSetup.jl")
using .BasiliskSetup

# Vizard interface (can now use BasiliskSetup)
include("VizardInterface.jl")
using .VizardInterface

# Export VisualizationCore
export SpacecraftGeometry, LightingConfig, CameraConfig, BackgroundConfig
export VideoExportSettings, GameVisualizationConfig
export cubesat_1u, cubesat_3u, cubesat_6u
export default_camera, default_lighting
export hd_video, uhd_video

# Export BasiliskSetup
export check_basilisk, install_basilisk, get_basilisk_info, remove_basilisk_venv

# Export VizardInterface
export VizardScenario, write_vizard_file, visualize_with_vizard, run_vizard

end # module Visualization