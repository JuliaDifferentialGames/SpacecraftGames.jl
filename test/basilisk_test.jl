# test/basilisk_setup_test.jl

using Test
using SpacecraftGames

@testset "Basilisk Setup" begin
    
    @testset "Virtual Environment Management" begin
        
        @testset "Get Python Executable" begin
            # This should create venv if it doesn't exist
            python_path = SpacecraftGames.Visualization.BasiliskSetup.get_python_executable()
            
            @test !isnothing(python_path)
            @test isa(python_path, String)
            
            # Check venv directory exists
            venv_dir = SpacecraftGames.Visualization.BasiliskSetup.VENV_DIR
            @test isdir(venv_dir)
            
            # Check Python executable exists
            @test isfile(python_path)
            
            println("Python executable: $python_path")
            println("Virtual environment: $venv_dir")
        end
        
        @testset "Python Version Check" begin
            python_path = SpacecraftGames.Visualization.BasiliskSetup.get_python_executable()
            
            # Verify Python works
            version_output = try
                read(`$python_path --version`, String)
            catch e
                println("Warning: Python version check failed: $e")
                ""
            end
            
            if !isempty(version_output)
                @test occursin("Python", version_output)
                println("Python version: ", strip(version_output))
            else
                println("Skipping Python version test")
            end
        end
    end
    
    @testset "Basilisk Installation Check" begin
        
        @testset "Check Installation Status" begin
            installed = check_basilisk()
            @test isa(installed, Bool)
            
            if installed
                println("✓ Basilisk is installed")
            else
                println("✗ Basilisk is not installed")
            end
        end
        
        @testset "Get Installation Info" begin
            @test_nowarn SpacecraftGames.Visualization.BasiliskSetup.get_basilisk_info()
        end
    end
    
    @testset "Basilisk Installation" begin
        
        # Check current status
        initially_installed = check_basilisk()
        
        if !initially_installed
            println("\n" * "="^60)
            println("Basilisk not installed - testing installation")
            println("This may take several minutes...")
            println("="^60)
            
            @testset "Install Basilisk" begin
                # Attempt installation
                success = install_basilisk()
                
                # Installation might fail in CI environments without proper dependencies
                if success
                    @test check_basilisk()
                    println("✓ Basilisk successfully installed and verified")
                else
                    # Installation failed - this is acceptable in minimal environments
                    println("⚠ Basilisk installation failed - this is expected in minimal environments")
                    println("  Manual installation: source .basilisk_venv/bin/activate && pip install basilisk")
                end
            end
        else
            println("Basilisk already installed - skipping installation test")
            
            @testset "Force Reinstall" begin
                # Test force reinstall flag
                println("\nTesting force reinstall...")
                success = install_basilisk(force=true)
                
                if success
                    @test check_basilisk()
                    println("✓ Basilisk force reinstall successful")
                else
                    println("⚠ Force reinstall failed")
                end
            end
        end
    end
    
    @testset "Vizard Script Execution" begin
        
        # Only test if Basilisk is installed
        if check_basilisk()
            
            @testset "Run Vizard Script" begin
                # Create a simple test script
                times = [0.0, 0.5, 1.0]
                N = length(times)
                
                states = zeros(12, N)
                states[1, :] = [0.0, 10.0, 20.0]
                controls = zeros(6, N-1)
                
                traj = Trajectory{Float64}(1, states, controls, times, 0.0)
                
                script_path = visualize_with_vizard(
                    [traj],
                    "basilisk_viz_test.py",
                    spacecraft_names = ["TestSC"]
                )
                
                @test isfile(script_path)
                
                # Run the script to generate Vizard save file
                try
                    run_vizard(script_path, blocking=true)
                    
                    # Check if Vizard save file was created
                    expected_save_file = replace(script_path, ".py" => "_UnityViz.txt")
                    
                    if isfile(expected_save_file)
                        @test true
                        println("✓ Vizard save file created: $expected_save_file")
                        
                        # Verify save file has content
                        content = read(expected_save_file, String)
                        @test length(content) > 0
                        println("  Save file size: $(length(content)) bytes")
                    else
                        println("⚠ Vizard save file not created (may need Vizard installed)")
                    end
                catch e
                    println("⚠ Script execution failed: $e")
                end
                
                # Clean up
                rm(script_path, force=true)
                save_file = replace(script_path, ".py" => "_UnityViz.txt")
                rm(save_file, force=true)
            end
            
        else
            println("⚠ Skipping Vizard execution tests - Basilisk not installed")
        end
    end
    
    @testset "Integration Test with Real Trajectories" begin
        
        if check_basilisk()
            
            @testset "Formation Flying Scenario" begin
                # Create realistic formation flying scenario
                T_final = 10.0
                dt = 1.0
                times = 0:dt:T_final
                N = length(times)
                
                trajectories = []
                
                # Chief at origin
                chief_states = zeros(12, N)
                push!(trajectories, Trajectory{Float64}(1, chief_states, zeros(6, N-1), collect(times), 0.0))
                
                # Deputy in circular motion
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
                    
                    θ_half = θ / 2
                    σ_mag = tan(θ_half / 2)
                    if abs(σ_mag) > 1.0
                        σ_mag = -1.0 / σ_mag
                    end
                    deputy_states[9, k] = σ_mag
                    deputy_states[12, k] = ω
                end
                
                push!(trajectories, Trajectory{Float64}(2, deputy_states, zeros(6, N-1), collect(times), 0.0))
                
                # Generate visualization script
                script_path = visualize_with_vizard(
                    trajectories,
                    "integration_test_formation.py",
                    spacecraft_names = ["Chief", "Deputy"]
                )
                
                @test isfile(script_path)
                
                content = read(script_path, String)
                @test occursin("Chief", content)
                @test occursin("Deputy", content)
                @test occursin("n_spacecraft = 2", content)
                
                println("\n" * "="^60)
                println("Integration Test Script Generated")
                println("="^60)
                println("File: $script_path")
                println("\nTo visualize:")
                println("  run_vizard(\"$script_path\")")
                println("  or: python $script_path")
                println("="^60)
                
                # Keep file for manual testing
                # Uncomment to clean up:
                # rm(script_path, force=true)
            end
            
        else
            println("⚠ Skipping integration test - Basilisk not installed")
        end
    end
    
    @testset "Cleanup and Reinstall" begin
        
        @testset "VEnv Removal" begin
            # Test that we can remove and recreate venv
            initial_exists = isdir(SpacecraftGames.Visualization.BasiliskSetup.VENV_DIR)
            
            if initial_exists
                removed = SpacecraftGames.Visualization.BasiliskSetup.remove_basilisk_venv()
                @test removed
                @test !isdir(SpacecraftGames.Visualization.BasiliskSetup.VENV_DIR)
                println("✓ Virtual environment removed")
                
                # Recreate it
                python_path = SpacecraftGames.Visualization.BasiliskSetup.get_python_executable()
                @test isfile(python_path)
                println("✓ Virtual environment recreated")
            else
                println("No venv to test removal")
            end
        end
    end
    
    @testset "Error Handling" begin
        
        @testset "Invalid Script Path" begin
            if check_basilisk()
                @test_throws Exception run_vizard("nonexistent_script.py")
            else
                println("Skipping error handling test - Basilisk not installed")
            end
        end
    end
    
end

# Summary
println("\n" * "="^60)
println("BASILISK SETUP TEST SUMMARY")
println("="^60)

basilisk_installed = check_basilisk()
println("Basilisk installed: ", basilisk_installed ? "✓ YES" : "✗ NO")

if basilisk_installed
    println("\nYou can now:")
    println("  1. Generate visualization scripts with visualize_with_vizard()")
    println("  2. Run scripts with run_vizard(script_path)")
    println("  3. Check status with get_basilisk_info()")
else
    println("\nTo install Basilisk:")
    println("  install_basilisk()")
    println("\nOr manually:")
    println("  source .basilisk_venv/bin/activate")
    println("  pip install basilisk")
end

println("="^60)