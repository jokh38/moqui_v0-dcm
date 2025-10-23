# MOQUI DCM Export Implementation Plan

## Review Summary

### TwoCentimeter Mode Analysis

**Configuration Location**: `moqui_tps.in` line 17
```
TwoCentimeterMode true
```

**Implementation Location**: `moqui/base/environments/mqi_tps_env.hpp`

**Key Implementation Details**:

1. **Variable Declaration** (line 108):
   ```cpp
   bool twoCentimeterMode = false;
   ```

2. **Configuration Parsing** (line 206):
   ```cpp
   this->twoCentimeterMode = parser.get_bool("twoCentimeterMode", false);
   ```

3. **Special Geometry Setup** (lines 787-800):
   - When `twoCentimeterMode` is true AND `usingPhantomGeo` is true:
   - Sets phantom dimensions to 400×1×400 (X×Y×Z)
   - Sets phantom units: X=1.0f, Y=2.0f, Z=1.0f (where Y units are in cm, X/Z in mm)
   - **Single voxel in Y direction = 2cm physical thickness** (overall phantom thickness)
   - Positions phantom centered at Y=-1.0
   - Creates **three separate phantom regions** in Z direction:
     - **Front phantom**: Z from 1mm to 20mm (19mm thick, 19 voxels)
     - **Center phantom**: Z from -1mm to 1mm (2 voxels, the measurement plane)
     - **Back phantom**: Z from -380mm to -1mm (189 voxels, 379mm thick)

4. **Purpose**: Creates a water slab phantom (400×400×2cm) for patient-specific QA dose measurements. The phantom slab lies in the XY plane with 2cm thickness in Y direction, but dose measurement plane is 2 voxels thick in Z direction at Z=0.

### Current Raw File Output Implementation

**Configuration**: `moqui_tps.in` line ~50
```
OutputFormat raw
```

**Implementation Files**:
- `moqui/base/mqi_io.hpp` - Contains actual file writing functions
- `moqui/base/environments/mqi_tps_env.hpp` - Contains output logic

**Key Output Functions**:
- `save_to_bin<R>()` - Main binary output function
- `save_reshaped_files()` - Handles reshaped data output
- `save_sparse_file()` - Handles NPZ format output

**Current Output Formats Supported**:
- `raw` - Binary format (default)
- `mhd` - MetaImage format
- `mha` - MetaImage format
- `npz` - NumPy compressed format

### Files That Need Revision for DCM Export

1. **Primary Files**:
   - `moqui/base/mqi_io.hpp` - Add DCM output functions
   - `moqui/base/environments/mqi_tps_env.hpp` - Add DCM format handling logic

2. **Configuration Files**:
   - `moqui_tps.in` - Already has `OutputFormat` parameter

3. **Build System**:
   - Need to add DCMTK library dependency to CMake configuration

### DCMTK Integration Requirements

**Dependencies to Add**:
- DCMTK libraries for DICOM file writing
- DICOM RT Dose format support

**Key Functions to Implement**:
- `save_to_dcm<R>()` - Main DICOM RT Dose export function
- DICOM header generation functions
- DICOM metadata handling (patient info, study info, dose scaling)

**Implementation Steps**:
1. Add DCMTK dependency to build system
2. Create DICOM RT Dose writer class
3. Implement dose data conversion to DICOM format
4. Add necessary DICOM metadata handling
5. Integrate with existing output pipeline
6. Test with TwoCentimeterMode output

### Integration Points

**TwoCentimeterMode Integration**:
- DCM export should work with both regular phantom and TwoCentimeterMode geometries
- Need to handle proper spatial scaling for 2cm slice
- Ensure DICOM coordinate system matches phantom geometry

**Output Pipeline Integration**:
- Add `OutputFormat dcm` option in `moqui_tps.in`
- Modify `save_reshaped_files()` to call DCM export when format is "dcm"
- Ensure proper dose scaling and units conversion