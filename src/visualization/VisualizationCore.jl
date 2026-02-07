# src/visualization/VisualizationCore.jl

"""
Core visualization types and interfaces for SpacecraftGames.jl
"""
module VisualizationCore

using LinearAlgebra
using StaticArrays
using Rotations
using ColorTypes

# Export types
export SpacecraftGeometry, LightingConfig, CameraConfig, BackgroundConfig
export VideoExportSettings, GameVisualizationConfig

# Export convenience constructors
export cubesat_1u, cubesat_3u, cubesat_6u
export default_camera, default_lighting
export hd_video, uhd_video

# Export utilities
export extract_position, extract_attitude, has_attitude, mrp_to_quaternion

# ============================================================================
# Spacecraft Geometry Configuration
# ============================================================================

"""
    SpacecraftGeometry

Physical and visual properties of a spacecraft for rendering.

# Fields
- `physical_dimensions::SVector{3,Float64}`: [length, width, height] in meters
- `visual_scale::Float64`: Rendering scale factor (>1.0 enlarges for visibility)
- `color::RGB{Float64}`: Spacecraft color
- `transparency::Float64`: Alpha value ∈ [0,1], 0=opaque, 1=transparent
- `label::String`: Optional label for identification

# Notes
For visibility at formation flying distances, visual_scale typically 10-50×.
Label in figure captions: "Spacecraft dimensions scaled {visual_scale}× for visibility"
"""
struct SpacecraftGeometry
    physical_dimensions::SVector{3,Float64}
    visual_scale::Float64
    color::RGB{Float64}
    transparency::Float64
    label::String
    
    function SpacecraftGeometry(
        dims::SVector{3,Float64},
        scale::Float64,
        color::RGB{Float64},
        transparency::Float64 = 0.0,
        label::String = ""
    )
        @assert all(dims .> 0) "Dimensions must be positive"
        @assert scale > 0 "Visual scale must be positive"
        @assert 0 ≤ transparency ≤ 1 "Transparency must be in [0,1]"
        new(dims, scale, color, transparency, label)
    end
end

"""
    cubesat_1u(; scale=15.0, color=RGB(0.8,0.8,0.9), label="")

Standard 1U CubeSat geometry (10cm cube).
"""
function cubesat_1u(; scale=15.0, color=RGB(0.8,0.8,0.9), transparency=0.0, label="")
    SpacecraftGeometry(SVector(0.1, 0.1, 0.1), scale, color, transparency, label)
end

"""
    cubesat_3u(; scale=15.0, color=RGB(0.8,0.8,0.9), label="")

Standard 3U CubeSat geometry (10cm × 10cm × 30cm).
"""
function cubesat_3u(; scale=15.0, color=RGB(0.8,0.8,0.9), transparency=0.0, label="")
    SpacecraftGeometry(SVector(0.1, 0.1, 0.3), scale, color, transparency, label)
end

"""
    cubesat_6u(; scale=15.0, color=RGB(0.8,0.8,0.9), label="")

Standard 6U CubeSat geometry (20cm × 10cm × 30cm).
"""
function cubesat_6u(; scale=15.0, color=RGB(0.8,0.8,0.9), transparency=0.0, label="")
    SpacecraftGeometry(SVector(0.2, 0.1, 0.3), scale, color, transparency, label)
end

# ============================================================================
# Lighting Configuration
# ============================================================================

"""
    LightingConfig

Space lighting model configuration.

# Fields
- `sun_direction_lvlh::SVector{3,Float64}`: Sun direction unit vector in LVLH frame
- `sun_color::RGB{Float64}`: Solar spectrum color (slightly warm)
- `ambient_intensity::Float64`: Non-physical fill light ∈ [0,1]
- `ambient_color::RGB{Float64}`: Fill light color (slight blue tint for space)
- `update_sun_direction::Bool`: Whether to animate sun direction over orbit

# Notes
Standard aerospace visualization uses ambient_intensity ≈ 0.2 for shadow fill.
This is non-physical but improves visibility. Document in captions:
"Lighting simplified for clarity - ambient illumination added"
"""
struct LightingConfig
    sun_direction_lvlh::SVector{3,Float64}
    sun_color::RGB{Float64}
    ambient_intensity::Float64
    ambient_color::RGB{Float64}
    update_sun_direction::Bool
    
    function LightingConfig(
        sun_dir::SVector{3,Float64},
        sun_color::RGB{Float64} = RGB(1.0, 0.95, 0.9),
        ambient_intensity::Float64 = 0.2,
        ambient_color::RGB{Float64} = RGB(0.2, 0.2, 0.25),
        update_sun::Bool = false
    )
        @assert isapprox(norm(sun_dir), 1.0, atol=1e-6) "Sun direction must be unit vector"
        @assert 0 ≤ ambient_intensity ≤ 1 "Ambient intensity must be in [0,1]"
        new(sun_dir, sun_color, ambient_intensity, ambient_color, update_sun)
    end
end

"""
    default_lighting(; sun_direction=SVector(1.0,0.0,0.0))

Default space lighting configuration with sun along +x LVLH.
"""
function default_lighting(; sun_direction=SVector(1.0, 0.0, 0.0))
    LightingConfig(normalize(sun_direction))
end

# ============================================================================
# Camera Configuration
# ============================================================================

"""
    CameraConfig

LVLH-locked camera configuration.

# Fields
- `offset_lvlh::SVector{3,Float64}`: Camera position offset from chief in LVLH [m]
- `up_direction::SVector{3,Float64}`: Camera "up" vector (typically -z = nadir)
- `fov::Float64`: Field of view in degrees

# Notes
Camera maintains fixed attitude relative to LVLH axes (Interpretation A).
Position tracks chief spacecraft, resulting in rotating view in inertial space.
Typical offset: behind and above chief, e.g., [-100, 0, 50]
"""
struct CameraConfig
    offset_lvlh::SVector{3,Float64}
    up_direction::SVector{3,Float64}
    fov::Float64
    
    function CameraConfig(
        offset::SVector{3,Float64},
        up::SVector{3,Float64} = SVector(0.0, 0.0, -1.0),
        fov::Float64 = 45.0
    )
        @assert isapprox(norm(up), 1.0, atol=1e-6) "Up direction must be unit vector"
        @assert 10 ≤ fov ≤ 120 "FOV should be reasonable (10-120 degrees)"
        new(offset, up, fov)
    end
end

"""
    default_camera(; distance=100.0, elevation=30.0)

Default camera positioned behind and above chief.

# Arguments
- `distance::Float64`: Distance behind chief along -x LVLH [m]
- `elevation::Float64`: Height above chief along +z LVLH [m]
"""
function default_camera(; distance=100.0, elevation=30.0)
    CameraConfig(SVector(-distance, 0.0, elevation))
end

# ============================================================================
# Background Configuration
# ============================================================================

"""
    BackgroundConfig

Starfield background configuration.

# Fields
- `show_stars::Bool`: Whether to render star background
- `star_density::Int`: Number of stars to generate
- `star_brightness::Float64`: Base star brightness ∈ [0,1]
- `background_color::RGB{Float64}`: Space background color (near-black)

# Notes
Stars are rendered as stationary in LVLH frame (rotate with orbit).
For physically accurate stationary stars, would need counter-rotation
at orbit rate - omitted for simplicity, acceptable for presentations.
"""
struct BackgroundConfig
    show_stars::Bool
    star_density::Int
    star_brightness::Float64
    background_color::RGB{Float64}
    
    function BackgroundConfig(
        show_stars::Bool = true,
        density::Int = 1000,
        brightness::Float64 = 0.8,
        bg_color::RGB{Float64} = RGB(0.01, 0.01, 0.02)
    )
        @assert density > 0 "Star density must be positive"
        @assert 0 ≤ brightness ≤ 1 "Brightness must be in [0,1]"
        new(show_stars, density, brightness, bg_color)
    end
end

# ============================================================================
# Video Export Settings
# ============================================================================

"""
    VideoExportSettings

Configuration for video file generation.

# Fields
- `resolution::Tuple{Int,Int}`: (width, height) in pixels
- `framerate::Int`: Frames per second
- `format::String`: Container format ("mp4", "mkv", "avi")
- `codec::String`: Video codec ("h264", "h265")
- `quality::Int`: Quality parameter (lower = better, 18-28 typical)

# Notes
Standard conference presentation settings:
- Resolution: (1920, 1080) for HD, (3840, 2160) for 4K
- Framerate: 30 FPS (cinema standard) or 60 FPS (smoother)
- Quality: 20-23 for good balance of file size and quality
"""
struct VideoExportSettings
    resolution::Tuple{Int,Int}
    framerate::Int
    format::String
    codec::String
    quality::Int
    
    function VideoExportSettings(
        resolution::Tuple{Int,Int} = (1920, 1080),
        framerate::Int = 30,
        format::String = "mp4",
        codec::String = "h264",
        quality::Int = 23
    )
        @assert all(resolution .> 0) "Resolution must be positive"
        @assert framerate > 0 "Framerate must be positive"
        @assert format in ["mp4", "mkv", "avi", "webm"] "Unsupported format"
        @assert 0 ≤ quality ≤ 51 "Quality must be in [0,51]"
        new(resolution, framerate, format, codec, quality)
    end
end

"""
    hd_video(; framerate=30)

HD (1920×1080) video export settings.
"""
hd_video(; framerate=30) = VideoExportSettings((1920, 1080), framerate)

"""
    uhd_video(; framerate=30)

4K UHD (3840×2160) video export settings.
"""
uhd_video(; framerate=30) = VideoExportSettings((3840, 2160), framerate)

# ============================================================================
# Complete Visualization Configuration
# ============================================================================

"""
    GameVisualizationConfig{T}

Complete configuration for differential game visualization.

# Fields
- `chief_id::Int`: Player ID of chief spacecraft (camera reference)
- `spacecraft_geometries::Vector{SpacecraftGeometry}`: Geometry for each player
- `camera::CameraConfig`: Camera positioning and orientation
- `lighting::LightingConfig`: Lighting model
- `background::BackgroundConfig`: Star field configuration
- `export_settings::VideoExportSettings`: Video output parameters

# Constructor
    GameVisualizationConfig(n_players::Int; kwargs...)

Creates configuration for n_players with sensible defaults.

# Keyword Arguments
- `chief_id::Int = 1`: Chief spacecraft player ID
- `geometries::Vector{SpacecraftGeometry}`: Custom geometries (defaults to 1U CubeSats)
- `camera::CameraConfig`: Custom camera (defaults to 100m behind, 30m above)
- `lighting::LightingConfig`: Custom lighting (defaults to sun along +x)
- `background::BackgroundConfig`: Custom background (defaults to stars enabled)
- `export_settings::VideoExportSettings`: Custom export (defaults to HD 30fps)
"""
struct GameVisualizationConfig{T}
    chief_id::Int
    spacecraft_geometries::Vector{SpacecraftGeometry}
    camera::CameraConfig
    lighting::LightingConfig
    background::BackgroundConfig
    export_settings::VideoExportSettings
    
    function GameVisualizationConfig{T}(
        chief_id::Int,
        geometries::Vector{SpacecraftGeometry},
        camera::CameraConfig,
        lighting::LightingConfig,
        background::BackgroundConfig,
        exporter::VideoExportSettings
    ) where T
        @assert chief_id > 0 "Chief ID must be positive"
        @assert length(geometries) > 0 "Must have at least one spacecraft geometry"
        new{T}(chief_id, geometries, camera, lighting, background, exporter)
    end
end

"""
    GameVisualizationConfig(n_players::Int; kwargs...)

Create visualization configuration with defaults.
"""
function GameVisualizationConfig{T}(
    n_players::Int;
    chief_id::Int = 1,
    geometries::Union{Vector{SpacecraftGeometry}, Nothing} = nothing,
    camera::CameraConfig = default_camera(),
    lighting::LightingConfig = default_lighting(),
    background::BackgroundConfig = BackgroundConfig(),
    export_settings::VideoExportSettings = hd_video()
) where T
    
    # Default geometries: colored 1U CubeSats
    if isnothing(geometries)
        colors = [
            RGB(0.3, 0.6, 0.9),  # Blue (chief)
            RGB(0.9, 0.3, 0.3),  # Red
            RGB(0.3, 0.9, 0.3),  # Green
            RGB(0.9, 0.9, 0.3),  # Yellow
            RGB(0.9, 0.3, 0.9),  # Magenta
        ]
        geometries = [
            cubesat_1u(color=colors[mod1(i, length(colors))], 
                      label="Player $i")
            for i in 1:n_players
        ]
    else
        @assert length(geometries) == n_players "Must provide geometry for each player"
    end
    
    GameVisualizationConfig{T}(chief_id, geometries, camera, lighting, 
                                background, export_settings)
end

# Convenience constructor without type parameter
GameVisualizationConfig(n_players::Int; kwargs...) = 
    GameVisualizationConfig{Float64}(n_players; kwargs...)

# ============================================================================
# State Extraction Utilities
# ============================================================================

"""
    extract_position(state::AbstractVector, has_attitude::Bool=true)

Extract position from state vector [x, y, z, vx, vy, vz, σ₁, σ₂, σ₃, ω₁, ω₂, ω₃].
Returns SVector{3} of [x, y, z] in LVLH frame.
"""
function extract_position(state::AbstractVector)
    SVector(state[1], state[2], state[3])
end

"""
    extract_attitude(state::AbstractVector)

Extract MRP attitude from state vector [x, y, z, vx, vy, vz, σ₁, σ₂, σ₃, ω₁, ω₂, ω₃].
Returns SVector{3} of MRP parameters.
"""
function extract_attitude(state::AbstractVector)
    @assert length(state) ≥ 9 "State must have attitude (length ≥ 9)"
    SVector(state[7], state[8], state[9])
end

"""
    has_attitude(state::AbstractVector)

Check if state vector includes attitude (length ≥ 9).
"""
has_attitude(state::AbstractVector) = length(state) ≥ 9

"""
    mrp_to_quaternion(σ::AbstractVector)

Convert Modified Rodrigues Parameters to quaternion.

# Arguments
- `σ::AbstractVector`: MRP vector [σ₁, σ₂, σ₃]

# Returns
- `QuatRotation`: Quaternion representation

# Notes
Implements Schaub & Junkins (2014) Eq. 3.141.
Assumes MRP shadow set switching has been applied (‖σ‖ ≤ 1).
"""
function mrp_to_quaternion(σ::AbstractVector)
    @assert length(σ) == 3 "MRP must be 3-dimensional"
    
    σ_norm_sq = sum(abs2, σ)
    
    # Warn if near singularity (should use shadow set)
    if σ_norm_sq > 1.0
        @warn "MRP norm $(√σ_norm_sq) > 1, consider shadow set switching" maxlog=1
    end
    
    denom = 1.0 + σ_norm_sq
    
    QuatRotation(
        (1.0 - σ_norm_sq) / denom,  # Scalar part (w)
        2.0 * σ[1] / denom,         # i
        2.0 * σ[2] / denom,         # j
        2.0 * σ[3] / denom          # k
    )
end

end