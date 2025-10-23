# DICOM RT Dose Export Implementation Plan - MOQUI TPS

## Overview

This document describes the implementation plan for adding DICOM RT Dose export capability to MOQUI. The implementation will integrate with the existing output pipeline to support `OutputFormat dcm` alongside the current formats (raw, mhd, mha, npz).

## Current Output Implementation Analysis

### Configuration

**Location**: `tps_env/moqui_tps.in` line 42
```
OutputFormat raw
```

**Supported values**: `raw`, `mhd`, `mha`, `npz`

### Implementation Files

**Primary implementation**: `moqui/base/mqi_io.hpp` lines 25-597

**Output orchestration**: `moqui/base/environments/mqi_tps_env.hpp`
- Lines 92-93: Format variable declaration
- Lines 279-287: Configuration parsing
- Lines 1583-1622: `save_reshaped_files()` function
- Lines 1626-1649: `save_sparse_file()` function
- Lines 1148-1153: Output flow control

## Existing Output Pipeline Architecture

### Format Selection Logic (mqi_tps_env.hpp:279-287)

```cpp
output_format = parser.get_string("OutputFormat", "raw");

if (strcasecmp(output_format.c_str(), "npz") == 0) {
    this->reshape_output = false;  // Keep sparse structure
    this->sparse_output  = true;
} else {
    this->reshape_output = true;   // Convert to dense 3D array
    this->sparse_output  = false;
}
```

**Key insight**: NPZ uses sparse data directly, while raw/mhd/mha first reshape to dense arrays.

### Output Flow Control (mqi_tps_env.hpp:1148-1153)

```cpp
if (this->reshape_output) {
    this->save_reshaped_files();      // Calls mqi_io.hpp functions
} else if (this->sparse_output) {
    this->save_sparse_file();         // NPZ-specific path
}
```

### Dense Format Output: save_reshaped_files() (mqi_tps_env.hpp:1583-1622)

**Step 1**: Convert sparse hash table to dense 3D array (lines 1561-1579)
```cpp
double* reshaped_data = new double[vol_size]();  // Allocate dense array
for (int ind = 0; ind < max_capacity; ind++) {
    if (data_[ind].key1 != empty && data_[ind].key2 != empty) {
        reshaped_data[data_[ind].key1] += data_[ind].value;  // Accumulate
    }
}
```

**Step 2**: Format-specific output (lines 1597-1617)
```cpp
if (!this->output_format.compare("mhd")) {
    mqi::io::save_to_mhd<R>(this->world->children[c_ind],
                            reshaped_data,
                            this->particles_per_history,  // ★ Scale factor
                            this->output_path,
                            filename,
                            vol_size);
} else if (!this->output_format.compare("mha")) {
    mqi::io::save_to_mha<R>(/* same parameters */);
} else {  // Default: raw format
    mqi::io::save_to_bin<double>(reshaped_data,
                                 this->particles_per_history,
                                 this->output_path,
                                 filename,
                                 vol_size);
}
```

**Key parameters available**:
- `this->world->children[c_ind]` - Phantom node (contains grid geometry)
- `reshaped_data` - Dense 1D array of dose values
- `this->particles_per_history` - Scale factor (line 36 in config)
- `this->output_path` - Output directory path
- `filename` - Beam name + child index + scorer name
- `vol_size` - Total number of voxels

## Existing Output Format Implementations

### 1. Binary Raw Format (mqi_io.hpp:162-200)

**Function**: `save_to_bin<double>()` - Dense array version

**Output files**: `{filename}.raw`

**What it writes** (lines 196-198):
```cpp
std::valarray<double> dest(src, length);
dest *= scale;  // Apply particles_per_history scaling
fid.write(reinterpret_cast<const char*>(&dest[0]), length * sizeof(double));
```

**No metadata**: Just raw binary doubles, no coordinate information.

### 2. MetaImage MHD Format (mqi_io.hpp:495-549)

**Function**: `save_to_mhd<R>()`

**Output files**:
- `{filename}.mhd` - Text header with metadata
- `{filename}.raw` - Binary data (double precision)

**Metadata written to .mhd header** (lines 518-537):
```cpp
ObjectType = Image
NDims = 3
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 0 1 0 0 0 1
Offset = {x0} {y0} {z0}           // Voxel center in mm
DimSize = {nx} {ny} {nz}          // Grid dimensions
ElementType = MET_DOUBLE
ElementSpacing = {dx} {dy} {dz}   // Voxel spacing in mm
ElementDataFile = {filename}.raw
```

**Grid information extraction** (lines 499-515):
```cpp
const mqi::grid3d<R>& grid = child->geo_;
int    nx = grid.dim_[0];
int    ny = grid.dim_[1];
int    nz = grid.dim_[2];
double dx = grid.voxel_size_[0];
double dy = grid.voxel_size_[1];
double dz = grid.voxel_size_[2];
double x0 = grid.x_[grid.dim_[0] / 2];  // Voxel center coordinates
double y0 = grid.y_[grid.dim_[1] / 2];
double z0 = grid.z_[grid.dim_[2] / 2];
```

### 3. MetaImage MHA Format (mqi_io.hpp:553-597)

**Function**: `save_to_mha<R>()`

**Output files**: `{filename}.mha` - Single file with header + data

**Metadata** (lines 576-593): Similar to MHD, but:
- Line 590: `HeaderSize = -1` (automatic calculation)
- Line 593: `ElementDataFile = LOCAL` (data embedded)

**Data embedding** (line 594):
```cpp
fid_header.write(reinterpret_cast<const char*>(&dest[0]), length * sizeof(double));
```

### 4. NPZ Sparse Format (mqi_io.hpp:251-491)

**Function**: `save_to_npz<R>()` - Three variants

**Output files**: `{filename}.npz` - ZIP archive containing:
- `indices.npy` - Column indices (CSR sparse format)
- `indptr.npy` - Row index pointers
- `shape.npy` - Matrix dimensions `[num_spots, vol_size]`
- `data.npy` - Non-zero dose values
- `format.npy` - String "csr"

**Key difference**: Preserves sparse structure, doesn't create dense array.

## Required Grid Geometry Information

### Grid3D Structure (accessed in mqi_io.hpp:499-515)

The `mqi::grid3d<R>` object provides:

| Property | Type | Example | Usage |
|----------|------|---------|-------|
| `dim_[0]` | int | 400 | X dimension (number of voxels) |
| `dim_[1]` | int | 400 | Y dimension |
| `dim_[2]` | int | 1 | Z dimension |
| `voxel_size_[0]` | double | 1.0 | X voxel spacing (mm) |
| `voxel_size_[1]` | double | 1.0 | Y voxel spacing (mm) |
| `voxel_size_[2]` | double | 2.0 | Z voxel spacing (mm) |
| `x_[i]` | double | -200 to +200 | X coordinates of voxel centers (mm) |
| `y_[i]` | double | -200 to +200 | Y coordinates |
| `z_[i]` | double | -1 to +1 | Z coordinates |

**Example for 2cm mode** (see [2cmmode.md](2cmmode.md)):
- Grid: 400×400×1 voxels
- Voxel size: 1×1×2 mm³
- Physical extent: 400×400×2 mm³
- Origin: Center at isocenter (0, 0, 0)

## DICOM RT Dose Export Implementation Plan

### Step 1: Add DCMTK Dependency to Build System

**Files to modify**:
- `CMakeLists.txt` - Add DCMTK package finding and linking

**Required libraries**:
```cmake
find_package(DCMTK REQUIRED)
include_directories(${DCMTK_INCLUDE_DIRS})
target_link_libraries(moqui ${DCMTK_LIBRARIES})
```

**Specific DCMTK modules needed**:
- `dcmdata` - DICOM data structures
- `dcmrt` - DICOM RT objects (RT Dose, RT Plan, RT Structure Set)

### Step 2: Implement DICOM RT Dose Writer Function

**Location**: Add to `moqui/base/mqi_io.hpp` after line 597

**Function signature** (following existing pattern):
```cpp
template<typename R>
void save_to_dcm(const mqi::node<R>*  child,
                 const double*        src,
                 const R              scale,
                 const std::string&   filepath,
                 const std::string&   filename,
                 const size_t         length);
```

**Function declaration** (add after line 98):
```cpp
/// Save dose data as DICOM RT Dose file
template<typename R>
void save_to_dcm(const mqi::node<R>*  child,
                 const double*        src,
                 const R              scale,
                 const std::string&   filepath,
                 const std::string&   filename,
                 const size_t         length);
```

### Step 3: DICOM RT Dose Implementation Details

**Implementation structure** (add after line 597 in mqi_io.hpp):

```cpp
template<typename R>
void mqi::io::save_to_dcm(const mqi::node<R>*  child,
                         const double*        src,
                         const R              scale,
                         const std::string&   filepath,
                         const std::string&   filename,
                         const size_t         length)
{
    // Step 1: Extract grid geometry (similar to save_to_mhd lines 499-515)
    const mqi::grid3d<R>& grid = child->geo_;
    int    nx = grid.dim_[0];
    int    ny = grid.dim_[1];
    int    nz = grid.dim_[2];
    double dx = grid.voxel_size_[0];
    double dy = grid.voxel_size_[1];
    double dz = grid.voxel_size_[2];
    double x0 = grid.x_[grid.dim_[0] / 2];
    double y0 = grid.y_[grid.dim_[1] / 2];
    double z0 = grid.z_[grid.dim_[2] / 2];

    // Step 2: Apply scaling and find dose range
    std::valarray<double> dest(src, length);
    dest *= scale;
    double max_dose = dest.max();
    double min_dose = dest.min();

    // Step 3: Convert to 16-bit unsigned integers with dose scaling
    // DICOM uses: dose_value = pixel_value * DoseGridScaling
    double dose_grid_scaling = max_dose / 65535.0;  // Max 16-bit value
    std::vector<uint16_t> pixel_data(length);
    for (size_t i = 0; i < length; i++) {
        pixel_data[i] = static_cast<uint16_t>(dest[i] / dose_grid_scaling);
    }

    // Step 4: Create DICOM RT Dose object
    DRTDoseIOD rtdose;

    // Step 5: Set required DICOM attributes
    // (Details in next section)

    // Step 6: Write to file
    rtdose.write(filepath + "/" + filename + ".dcm");
}
```

### Step 4: Required DICOM RT Dose Attributes

**Module: Patient** (Type 2 - Required, may be empty)
```cpp
rtdose.getPatientName().putString("PHANTOM");
rtdose.getPatientID().putString("QA_PHANTOM");
rtdose.getPatientBirthDate().putString("");
rtdose.getPatientSex().putString("");
```

**Module: General Study** (Type 2)
```cpp
rtdose.getStudyInstanceUID().putString(generate_uid());
rtdose.getStudyDate().putString(get_current_date());
rtdose.getStudyTime().putString(get_current_time());
rtdose.getStudyID().putString("1");
```

**Module: RT Series** (Type 1 - Required)
```cpp
rtdose.getModality().putString("RTDOSE");
rtdose.getSeriesInstanceUID().putString(generate_uid());
rtdose.getSeriesNumber().putString("1");
```

**Module: Frame of Reference** (Type 1)
```cpp
rtdose.getFrameOfReferenceUID().putString(generate_uid());
rtdose.getPositionReferenceIndicator().putString("");
```

**Module: General Equipment** (Type 3 - Optional)
```cpp
rtdose.getManufacturer().putString("MOQUI");
rtdose.getSoftwareVersions().putString("1.0");
```

**Module: RT Dose** (Type 1)
```cpp
// Image dimensions
rtdose.getRows().putUint16(ny);
rtdose.getColumns().putUint16(nx);
rtdose.getNumberOfFrames().putString(std::to_string(nz).c_str());

// Pixel data properties
rtdose.getSamplesPerPixel().putUint16(1);
rtdose.getPhotometricInterpretation().putString("MONOCHROME2");
rtdose.getBitsAllocated().putUint16(16);
rtdose.getBitsStored().putUint16(16);
rtdose.getHighBit().putUint16(15);
rtdose.getPixelRepresentation().putUint16(0);  // Unsigned

// Dose attributes
rtdose.getDoseUnits().putString("GY");  // Gray
rtdose.getDoseType().putString("PHYSICAL");
rtdose.getDoseSummationType().putString("PLAN");
rtdose.getDoseGridScaling().putFloat64(dose_grid_scaling);

// Spatial attributes
OFString pixel_spacing;
pixel_spacing = std::to_string(dy) + "\\" + std::to_string(dx);
rtdose.getPixelSpacing().putString(pixel_spacing.c_str());

OFString grid_frame_offset_vector = "";
for (int i = 0; i < nz; i++) {
    if (i > 0) grid_frame_offset_vector += "\\";
    grid_frame_offset_vector += std::to_string(grid.z_[i] - z0);
}
rtdose.getGridFrameOffsetVector().putString(grid_frame_offset_vector.c_str());

// Image position (top-left corner of first voxel)
OFString image_position;
image_position = std::to_string(grid.x_[0] - dx/2) + "\\" +
                 std::to_string(grid.y_[0] - dy/2) + "\\" +
                 std::to_string(grid.z_[0] - dz/2);
rtdose.getImagePositionPatient().putString(image_position.c_str());

// Image orientation (default: HFS - Head First Supine)
rtdose.getImageOrientationPatient().putString("1\\0\\0\\0\\1\\0");

// Pixel data
rtdose.getPixelData().putUint16Array(pixel_data.data(), length);
```

### Step 5: Integration with save_reshaped_files()

**Location**: `moqui/base/environments/mqi_tps_env.hpp` lines 1597-1617

**Add DCM format handling**:
```cpp
if (!this->output_format.compare("mhd")) {
    mqi::io::save_to_mhd<R>(this->world->children[c_ind],
                            reshaped_data,
                            this->particles_per_history,
                            this->output_path,
                            filename,
                            vol_size);
} else if (!this->output_format.compare("mha")) {
    mqi::io::save_to_mha<R>(this->world->children[c_ind],
                            reshaped_data,
                            this->particles_per_history,
                            this->output_path,
                            filename,
                            vol_size);
} else if (!this->output_format.compare("dcm")) {  // ★ ADD THIS BLOCK
    mqi::io::save_to_dcm<R>(this->world->children[c_ind],
                            reshaped_data,
                            this->particles_per_history,
                            this->output_path,
                            filename,
                            vol_size);
} else {
    mqi::io::save_to_bin<double>(reshaped_data,
                                 this->particles_per_history,
                                 this->output_path,
                                 filename,
                                 vol_size);
}
```

### Step 6: Configuration Update

**No changes needed** - `OutputFormat` parameter already exists in `moqui_tps.in` line 42.

**Usage**:
```
OutputFormat dcm
```

### Step 7: Helper Functions to Implement

**UID Generation** (DICOM requires unique identifiers):
```cpp
std::string generate_uid() {
    // Format: root.timestamp.random
    // Example: 1.2.826.0.1.3680043.10.511.3.12345678901234567890
    std::string root = "1.2.826.0.1.3680043.10.511";  // Registered OID
    // Append timestamp and random number
    return root + "." + get_timestamp() + "." + get_random();
}
```

**Date/Time Functions**:
```cpp
std::string get_current_date() {
    // Format: YYYYMMDD
    time_t now = time(0);
    struct tm tstruct = *localtime(&now);
    char buf[9];
    strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
    return buf;
}

std::string get_current_time() {
    // Format: HHMMSS
    time_t now = time(0);
    struct tm tstruct = *localtime(&now);
    char buf[7];
    strftime(buf, sizeof(buf), "%H%M%S", &tstruct);
    return buf;
}
```

## Testing Strategy

### Test Case 1: 2cm Mode Output

**Configuration**: Use existing [2cmmode.md](2cmmode.md) setup
```
TwoCentimeterMode true
OutputFormat dcm
```

**Expected output**:
- File: `{beam_name}_1_{scorer_name}.dcm`
- Dimensions: 400×400×1 (Rows=400, Columns=400, NumberOfFrames=1)
- Pixel spacing: 1.0\1.0 mm
- Grid frame offset: 0.0 (single slice at Z=0)
- Dose units: GY
- Data type: 16-bit unsigned integers

**Validation**:
1. Open in DICOM viewer (3D Slicer, MIM, RayStation)
2. Verify dimensions match expected 400×400
3. Verify dose values are scaled correctly
4. Verify spatial coordinates (ImagePositionPatient, PixelSpacing)

### Test Case 2: Regular Patient CT Output

**Configuration**: Full patient CT phantom
```
TwoCentimeterMode false
OutputFormat dcm
```

**Expected output**:
- Multiple .dcm files if multiple scorers exist
- Variable dimensions based on CT grid
- Proper coordinate transformation from CT space

**Validation**:
1. Load both CT and dose in treatment planning system
2. Verify dose overlays correctly on CT anatomy
3. Check coordinate system alignment

### Test Case 3: Compare with MHD Output

**Test**: Run same beam plan with both formats
```
OutputFormat mhd  # First run
OutputFormat dcm  # Second run
```

**Validation**:
1. Load MHD and DCM in same viewer
2. Subtract dose distributions
3. Verify difference is within numerical precision (< 0.001%)

## Implementation Checklist

- [ ] **CMake**: Add DCMTK dependency
- [ ] **mqi_io.hpp**: Add `save_to_dcm()` function declaration (after line 98)
- [ ] **mqi_io.hpp**: Implement `save_to_dcm()` function (after line 597)
- [ ] **mqi_io.hpp**: Add helper functions (UID generation, date/time)
- [ ] **mqi_tps_env.hpp**: Add DCM format handling in `save_reshaped_files()` (lines 1597-1617)
- [ ] **Build**: Compile with DCMTK, verify no linking errors
- [ ] **Test**: 2cm mode with DCM output
- [ ] **Test**: Regular phantom with DCM output
- [ ] **Test**: Compare DCM vs MHD dose distributions
- [ ] **Documentation**: Update user guide with DCM format usage

## File Summary

| File | Lines | Changes Required |
|------|-------|------------------|
| `CMakeLists.txt` | TBD | Add DCMTK package and linking |
| `moqui/base/mqi_io.hpp` | After 98 | Add function declaration |
| `moqui/base/mqi_io.hpp` | After 597 | Implement save_to_dcm() (~150 lines) |
| `moqui/base/mqi_io.hpp` | After save_to_dcm | Add helper functions (~50 lines) |
| `moqui/base/environments/mqi_tps_env.hpp` | 1597-1617 | Add else-if block for DCM (~10 lines) |
| `tps_env/moqui_tps.in` | Line 42 | No changes (parameter exists) |

## Related Documentation

- **2cm Mode Details**: See [2cmmode.md](2cmmode.md) for phantom geometry and scorer configuration
- **DICOM Standard**: DICOM PS3.3 Section C.8.8 (RT Dose Module)
- **DCMTK Documentation**: https://support.dcmtk.org/docs/
- **Current Output Formats**: Lines 25-597 in `moqui/base/mqi_io.hpp`