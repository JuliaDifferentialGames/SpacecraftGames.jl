# test/vizard_test.jl

using Test
using SpacecraftGames
using LinearAlgebra
using StaticArrays

@testset "Vizard Interface" begin
    
    @testset "Trajectory Conversion" begin
        # Create simple test trajectory
        T_final = 10.0
        dt = 1.0
        times = 0:dt:T_final
        N = length(times)
        
        # Chief spacecraft
        chief_states = zeros(12, N)
        chief_controls = zeros(6, N-1)
        
        traj_chief = Trajectory{Float64}(1, chief_states, chief_controls, collect(times), 0.0)
        
        # Deputy spacecraft with circular motion
        deputy_states = zeros(12, N)
        radius = 50.0
        ω = 2π / T_final
        
        for (k, t) in enumerate(times)
            θ = ω * t
            deputy_states[1, k] = radius * cos(θ)
            deputy_states[2, k] = radius * sin(θ)
            deputy_states[3, k] = 0.0
            deputy_states[4, k] = -radius * ω * sin(θ)
            deputy_states[5, k] = radius * ω * cos(θ)
            deputy_states[6, k] = 0.0
            
            # Attitude with shadow set
            θ_half = θ / 2
            σ_mag = tan(θ_half / 2)
            if abs(σ_mag) > 1.0
                σ_mag = -1.0 / σ_mag
            end
            deputy_states[9, k] = σ_mag
            deputy_states[12, k] = ω
        end
        
        traj_deputy = Trajectory{Float64}(2, deputy_states, zeros(6, N-1), collect(times), 0.0)
        
        trajectories = [traj_chief, traj_deputy]
        
        # Test conversion
        @testset "VizardScenario Creation" begin
            scenario = SpacecraftGames.Visualization.VizardInterface.trajectory_to_vizard(
                trajectories,
                ["Chief", "Deputy"]
            )
            
            @test length(scenario.spacecraft_data) == 2
            @test length(scenario.time_data) == N
            @test scenario.settings["n_spacecraft"] == 2
            @test scenario.settings["n_timesteps"] == N
            @test scenario.settings["reference_frame"] == "LVLH"
            
            # Check chief data
            chief_data = scenario.spacecraft_data[1]
            @test chief_data["name"] == "Chief"
            @test chief_data["player_id"] == 1
            @test size(chief_data["position"]) == (3, N)
            @test size(chief_data["velocity"]) == (3, N)
            @test size(chief_data["mrp"]) == (3, N)
            @test chief_data["has_attitude"] == true
            
            # Check deputy data
            deputy_data = scenario.spacecraft_data[2]
            @test deputy_data["name"] == "Deputy"
            @test deputy_data["player_id"] == 2
            
            # Verify position data is correct
            @test deputy_data["position"][1, 1] ≈ radius  # x at t=0
            @test deputy_data["position"][2, 1] ≈ 0.0     # y at t=0
        end
    end
    
    @testset "Python Script Generation" begin
        # Create minimal trajectory
        times = [0.0, 1.0, 2.0]
        N = length(times)
        
        states = zeros(12, N)
        states[1, :] = [0.0, 10.0, 20.0]  # x positions
        
        controls = zeros(6, N-1)
        traj = Trajectory{Float64}(1, states, controls, times, 0.0)
        
        scenario = SpacecraftGames.Visualization.VizardInterface.trajectory_to_vizard(
            [traj],
            ["TestSC"]
        )
        
        # Write to file
        output_file = "test_vizard_output.py"
        
        @testset "File Writing" begin
            result = SpacecraftGames.Visualization.VizardInterface.write_vizard_file(
                scenario,
                output_file
            )
            
            @test isfile(result)
            @test endswith(result, ".py")
            
            # Read file and check contents
            content = read(result, String)
            
            # Check for required imports
            @test occursin("import numpy as np", content)
            @test occursin("from Basilisk.utilities import vizSupport", content)
            
            # Check for scenario data
            @test occursin("n_spacecraft = 1", content)
            @test occursin("n_timesteps = 3", content)
            @test occursin("reference_frame = 'LVLH'", content)
            
            # Check for time data
            @test occursin("times = np.array([", content)
            @test occursin("0.0,", content)
            @test occursin("1.0,", content)
            @test occursin("2.0", content)
            
            # Check for spacecraft data
            @test occursin("# Spacecraft 1: TestSC", content)
            @test occursin("sc1_position", content)
            @test occursin("sc1_velocity", content)
            @test occursin("sc1_mrp", content)
            
            # Check for Vizard setup
            @test occursin("def run_visualization():", content)
            @test occursin("SimulationBaseClass.SimBaseClass()", content)
            @test occursin("vizSupport.enableUnityVisualization", content)
            
            # Clean up
            rm(result, force=true)
        end
    end
    
    @testset "High-level Interface" begin
        # Create test trajectories
        times = [0.0, 0.5, 1.0]
        N = length(times)
        
        trajectories = []
        
        for i in 1:3
            states = zeros(12, N)
            states[1, :] = [0.0, i*10.0, i*20.0]  # Different x positions
            controls = zeros(6, N-1)
            push!(trajectories, Trajectory{Float64}(i, states, controls, times, 0.0))
        end
        
        @testset "With Custom Names" begin
            output_file = "test_custom_names.py"
            
            result = visualize_with_vizard(
                trajectories,
                output_file,
                spacecraft_names = ["Alpha", "Beta", "Gamma"]
            )
            
            @test isfile(result)
            
            content = read(result, String)
            @test occursin("Alpha", content)
            @test occursin("Beta", content)
            @test occursin("Gamma", content)
            @test occursin("n_spacecraft = 3", content)
            
            rm(result, force=true)
        end
        
        @testset "With Default Names" begin
            output_file = "test_default_names.py"
            
            result = visualize_with_vizard(trajectories, output_file)
            
            @test isfile(result)
            
            content = read(result, String)
            @test occursin("SC1", content)
            @test occursin("SC2", content)
            @test occursin("SC3", content)
            
            rm(result, force=true)
        end
        
        @testset "Default Filename" begin
            result = visualize_with_vizard(trajectories)
            
            @test isfile(result)
            @test endswith(result, ".py")
            
            rm(result, force=true)
        end
    end
    
    @testset "Formation Flying Example" begin
        # Create a realistic formation flying scenario
        T_final = 20.0
        dt = 0.1
        times = 0:dt:T_final
        N = length(times)
        
        trajectories = []
        
        # Chief at origin
        chief_states = zeros(12, N)
        push!(trajectories, Trajectory{Float64}(1, chief_states, zeros(6, N-1), collect(times), 0.0))
        
        # Two deputies in formation
        for i in 1:2
            deputy_states = zeros(12, N)
            radius = 50.0
            phase = (i-1) * π  # 180° apart
            ω = 2π / T_final
            
            for (k, t) in enumerate(times)
                θ = ω * t + phase
                deputy_states[1, k] = radius * cos(θ)
                deputy_states[2, k] = radius * sin(θ)
                deputy_states[3, k] = 0.0
                deputy_states[4, k] = -radius * ω * sin(θ)
                deputy_states[5, k] = radius * ω * cos(θ)
                deputy_states[6, k] = 0.0
                
                θ_half = θ / 2
                σ_mag = tan(θ_half / 2)
                if abs(σ_mag) > 1.0
                    σ_mag = -1.0 / σ_mag
                end
                deputy_states[9, k] = σ_mag
                deputy_states[12, k] = ω
            end
            
            push!(trajectories, Trajectory{Float64}(i+1, deputy_states, zeros(6, N-1), collect(times), 0.0))
        end
        
        # Generate visualization
        output_file = "formation_example.py"
        result = visualize_with_vizard(
            trajectories,
            output_file,
            spacecraft_names = ["Chief", "Deputy1", "Deputy2"]
        )
        
        @test isfile(result)
        
        content = read(result, String)
        
        # Verify all spacecraft are included
        @test occursin("Chief", content)
        @test occursin("Deputy1", content)
        @test occursin("Deputy2", content)
        
        # Verify trajectory data structure
        @test occursin("sc1_position", content)
        @test occursin("sc2_position", content)
        @test occursin("sc3_position", content)
        
        # Verify time array is correct length
        @test occursin("n_timesteps = $(N)", content)
        
        println("\n" * "="^60)
        println("Formation flying example generated: $result")
        println("To visualize (requires Basilisk):")
        println("  python $result")
        println("="^60)
        
        # Keep this file for manual testing
        # Uncomment to clean up:
        # rm(result, force=true)
    end
    
end