# src/visualization/BasiliskSetup.jl

"""
Basilisk installation and management for SpacecraftGames.jl

Handles automatic installation of Basilisk in a Python virtual environment.
"""
module BasiliskSetup

# using Pkg

export check_basilisk, install_basilisk, get_python_executable, run_vizard_script

# Path to virtual environment within package
const VENV_DIR = joinpath(@__DIR__, "..", "..", ".basilisk_venv")
const VENV_PYTHON = Sys.iswindows() ? 
    joinpath(VENV_DIR, "Scripts", "python.exe") :
    joinpath(VENV_DIR, "bin", "python")

"""
    get_python_executable()

Get path to Python executable in virtual environment.
Creates venv if it doesn't exist.
"""
function get_python_executable()
    if !isdir(VENV_DIR)
        @info "Creating Python virtual environment at $VENV_DIR"
        create_venv()
    end
    return VENV_PYTHON
end

"""
    create_venv()

Create a Python virtual environment for Basilisk.
"""
function create_venv()
    # Check if python3 or python is available
    python_cmd = nothing
    for cmd in ["python3", "python"]
        try
            run(pipeline(`$cmd --version`, stdout=devnull, stderr=devnull))
            python_cmd = cmd
            break
        catch
            continue
        end
    end
    
    if python_cmd === nothing
        error("Python not found. Please install Python 3.7+ and ensure it's in your PATH.")
    end
    
    # Create virtual environment
    try
        run(`$python_cmd -m venv $VENV_DIR`)
        @info "✓ Virtual environment created at $VENV_DIR"
    catch e
        error("Failed to create virtual environment. Error: $e\n" *
              "You may need to install python3-venv: sudo apt-get install python3-venv")
    end
end

"""
    check_basilisk()

Check if Basilisk is installed in the virtual environment.

# Returns
- `true` if Basilisk is installed and importable
- `false` otherwise
"""
function check_basilisk()
    if !isfile(VENV_PYTHON)
        return false
    end
    
    # Try to import Basilisk (try multiple import patterns)
    test_script = """
import sys
try:
    # Try modern bsk package
    from Basilisk.utilities import SimulationBaseClass
    from Basilisk.utilities import vizSupport
    print("SUCCESS: Basilisk.utilities found")
    sys.exit(0)
except ImportError:
    try:
        # Try alternative import
        import Basilisk
        print("SUCCESS: Basilisk module found")
        sys.exit(0)
    except ImportError:
        try:
            # Try bsk package directly
            import bsk
            print("SUCCESS: bsk package found")
            sys.exit(0)
        except ImportError as e:
            print(f"IMPORT_ERROR: {e}")
            sys.exit(1)
"""
    
    try
        result = read(`$VENV_PYTHON -c $test_script`, String)
        return occursin("SUCCESS", result)
    catch e
        return false
    end
end

# In src/visualization/BasiliskSetup.jl, update install_basilisk():

function install_basilisk(; force::Bool=false)
    
    # Check if already installed
    if !force && check_basilisk()
        @info "Basilisk is already installed"
        return true
    end
    
    # Ensure venv exists
    python = get_python_executable()
    
    @info "Installing Basilisk..."
    println("This may take several minutes...")
    
    # Upgrade pip first
    try
        @info "Upgrading pip..."
        run(`$python -m pip install --upgrade pip setuptools wheel`)
    catch e
        @warn "Failed to upgrade pip: $e"
    end
    
    # Install Basilisk
    try
        @info "Installing basilisk package (bsk)..."
        # The correct package is 'bsk' on PyPI
        cmd = `$python -m pip install bsk`
        
        # Run with output visible
        process = run(cmd, wait=true)
        
        if process.exitcode == 0
            @info "✓ Basilisk installation complete"
            
            # Verify installation
            if check_basilisk()
                @info "✓ Basilisk verified and ready to use"
                return true
            else
                @warn "Package installed but verification failed"
                @warn "Checking what was installed..."
                
                # Debug: show what packages are installed
                try
                    pkg_list = read(`$python -m pip list`, String)
                    if occursin("bsk", pkg_list)
                        println("  bsk package is installed")
                    end
                catch
                end
                
                return false
            end
        else
            @error "Basilisk installation failed with exit code $(process.exitcode)"
            return false
        end
        
    catch e
        @error "Failed to install Basilisk" exception=e
        println("\nManual installation:")
        println("  1. Activate venv: source $VENV_DIR/bin/activate")
        println("  2. Install: pip install bsk")
        return false
    end
end

"""
    ensure_basilisk()

Ensure Basilisk is installed, installing if necessary.
Called automatically by visualization functions.
"""
function ensure_basilisk()
    if !check_basilisk()
        @info "Basilisk not found. Installing..."
        success = install_basilisk()
        if !success
            error("Failed to install Basilisk. Please install manually:\n" *
                  "  pip install basilisk")
        end
    end
end

"""
    run_vizard_script(script_path::String; blocking::Bool=true)

Run a Vizard Python script using the virtual environment.

# Arguments
- `script_path::String`: Path to the Python script
- `blocking::Bool`: If true, wait for script to complete

# Returns
- Process object if non-blocking
- Nothing if blocking
"""
function run_vizard_script(script_path::String; blocking::Bool=true)
    
    ensure_basilisk()
    
    if !isfile(script_path)
        error("Script not found: $script_path")
    end
    
    python = get_python_executable()
    
    @info "Running Vizard visualization..."
    @info "Script: $script_path"
    
    if blocking
        try
            run(`$python $script_path`)
            @info "✓ Vizard script completed"
            return nothing
        catch e
            @error "Vizard script failed" exception=e
            rethrow(e)
        end
    else
        # Run in background
        process = run(`$python $script_path`, wait=false)
        @info "Vizard running in background (PID: $(process.pid))"
        return process
    end
end

"""
    get_basilisk_info()

Get information about Basilisk installation.
"""
function get_basilisk_info()
    println("="^60)
    println("Basilisk Installation Information")
    println("="^60)
    
    println("Virtual environment: $VENV_DIR")
    println("Python executable: $VENV_PYTHON")
    println("VEnv exists: ", isdir(VENV_DIR))
    println("Python exists: ", isfile(VENV_PYTHON))
    
    if isfile(VENV_PYTHON)
        # Get Python version
        try
            version = read(`$VENV_PYTHON --version`, String)
            println("Python version: ", strip(version))
        catch
            println("Python version: Unable to determine")
        end
        
        # Check Basilisk
        installed = check_basilisk()
        println("Basilisk installed: ", installed)
        
        if installed
            # Get Basilisk version
            version_script = """
try:
    import Basilisk
    if hasattr(Basilisk, '__version__'):
        print(Basilisk.__version__)
    else:
        print("Version unknown")
except:
    print("Unable to determine")
"""
            try
                bsk_version = read(`$VENV_PYTHON -c $version_script`, String)
                println("Basilisk version: ", strip(bsk_version))
            catch
                println("Basilisk version: Unable to determine")
            end
        end
    else
        println("Basilisk installed: N/A (no Python)")
    end
    
    println("="^60)
end

"""
    remove_basilisk_venv()

Remove the Basilisk virtual environment.
Useful for clean reinstallation.
"""
function remove_basilisk_venv()
    if isdir(VENV_DIR)
        @info "Removing virtual environment at $VENV_DIR"
        rm(VENV_DIR, recursive=true, force=true)
        @info "✓ Virtual environment removed"
        return true
    else
        @info "No virtual environment to remove"
        return false
    end
end

end # module BasiliskSetup