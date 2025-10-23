# MOQUI v0 DCM Project Build Architecture

## Main Entry Point
- **File**: `tps_env/tps_env.cpp`
- **Purpose**: Main executable for running Monte Carlo simulations using the Treatment Planning System (TPS) environment
- **Functionality**: 
  - Initializes `mqi::tps_env` from input file
  - Runs the simulation
  - Reports total execution time
  - Command-line interface for the Moqui simulation engine

## Build Configuration
- **Main Executable**: `tps_env.cpp` links against the entire C++ codebase in `moqui/`
- **Configuration File**: `tps_env/moqui_tps.in` (default configuration)
- **Alternative Config**: Can be specified as command-line argument

## Configuration File Details (moqui_tps.in)
- GPU configuration (GPUID 0)
- Random seed setting
- Threading configuration (TotalThreads -1 for optimized)
- Batch size settings (MaxHistoriesPerBatch 10000000)
- Phantom geometry settings
- Dimensions: PhantomDimX 300, PhantomDimY 400, PhantomDimZ 300
- Scoring configuration (Dose)
- Input/output directories
- DICOM integration settings

## Key Dependencies
- GDCM library for DICOM handling
- CUDA support for GPU acceleration
- Moqui simulation engine headers in `moqui/base/`

## Build Process
The `tps_env.cpp` serves as the main entry point that links against all the Moqui simulation engine components located in the `moqui/` directory structure.