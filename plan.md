# DICOM RT Dose Export Implementation Validation - MOQUI TPS

## Overview

This document validates the **completed** DICOM RT Dose export implementation in MOQUI. The implementation uses **GDCM library** and is fully integrated with the existing output pipeline, supporting `OutputFormat dcm` alongside the current formats (raw, mhd, mha, npz).

**Status**: ✅ Implementation complete (see [mqi_io.hpp:637-856](moqui/base/mqi_io.hpp#L637-L856))
**Library**: GDCM (Grassroots DICOM)
**Integration**: Complete in [mqi_tps_env.hpp:1611-1617](moqui/base/environments/mqi_tps_env.hpp#L1611-L1617)

✅ **String Formatting Bugs Fixed**: All 4 critical string formatting bugs in separator characters (lines 805, 812, 822-824, 831) have been fixed.

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

## DICOM RT Dose Export Implementation Details

### Step 1: Build System Integration - ✅ COMPLETE

**Status**: GDCM library is already integrated in the build system.

**Current configuration** ([tps_env/CMakeLists.txt:37-68](tps_env/CMakeLists.txt#L37-L68)):
```cmake
find_package(GDCM REQUIRED)
include(${GDCM_USE_FILE})

include_directories(
    ${GDCM_INCLUDE_DIRS}
    # ...
)

target_link_libraries(tps_env PRIVATE
    gdcmCommon
    gdcmDSED
    gdcmMEXD
    gdcmjpeg12
    gdcmjpeg8
    gdcmDICT
    gdcmIOD
    gdcmMSFF
    gdcmjpeg16
    gdcmRT      # ← RT Dose support
)
```

**GDCM modules used**:
- `gdcmCommon` - Core DICOM data structures
- `gdcmDICT` - DICOM data dictionary
- `gdcmIOD` - Information Object Definitions
- `gdcmRT` - DICOM RT objects (RT Dose, RT Plan, RT Structure Set)

### Step 2: DICOM RT Dose Writer Function - ✅ COMPLETE

**Status**: Function is fully implemented and integrated.

**Declaration location**: [moqui/base/mqi_io.hpp:30-38](moqui/base/mqi_io.hpp#L30-L38) (in global `mqi::` namespace)

**Function signature** (actual implementation):
```cpp
template<typename R>
void
save_to_dcm(const mqi::node_t<R>* children,  // Note: node_t, not node
            const double*         src,
            const R               scale,
            const std::string&    filepath,
            const std::string&    filename,
            const uint32_t        length);    // Note: uint32_t, not size_t
```

**Implementation location**: [moqui/base/mqi_io.hpp:637-856](moqui/base/mqi_io.hpp#L637-L856)

**Note**: Uses `node_t<R>*` (not `node<R>*`) and parameter is named `children` (plural) to match existing API pattern.

### Step 3: DICOM RT Dose Implementation Details - ✅ COMPLETE

**Implementation structure** (actual code at [mqi_io.hpp:637-856](moqui/base/mqi_io.hpp#L637-L856)):

```cpp
template<typename R>
void mqi::io::save_to_dcm(const mqi::node_t<R>* children,  // ← Corrected type
                          const double*         src,
                          const R               scale,
                          const std::string&    filepath,
                          const std::string&    filename,
                          const uint32_t        length)      // ← Corrected type
{
    // Step 1: Extract grid geometry (similar to save_to_mhd)
    const mqi::grid3d<R>& grid = children->geo_[0];  // ← Array access [0]
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

    // Step 4: Create DICOM file using GDCM (not DCMTK)
    gdcm::File dcm_file;
    gdcm::DataSet& ds = dcm_file.GetDataSet();

    // Step 5: Set required DICOM attributes using GDCM Attribute API
    // (Details in next section - uses gdcm::Attribute<Group, Element>)

    // Step 6: Write to file using GDCM Writer
    gdcm::Writer writer;
    writer.SetFile(dcm_file);
    std::string full_path = filepath + "/" + filename + ".dcm";
    if (!writer.Write(full_path.c_str())) {
        std::cout << "Error: Failed to write DICOM RT Dose file" << std::endl;
    }
}
```

### Step 4: Required DICOM RT Dose Attributes - ✅ COMPLETE

**Implementation uses GDCM Attribute API** (not DCMTK):
- Pattern: `gdcm::Attribute<Group, Element> attr_name;`
- Setting: `attr_name.Set(value);`
- Inserting: `ds.Insert(attr_name.GetAsDataElement());`

**Module: Patient** (Type 2 - Required, may be empty)
```cpp
gdcm::Attribute<0x0010, 0x0010> patient_name;
patient_name.Set("PHANTOM");
ds.Insert(patient_name.GetAsDataElement());

gdcm::Attribute<0x0010, 0x0020> patient_id;
patient_id.Set("QA_PHANTOM");
ds.Insert(patient_id.GetAsDataElement());
```

**Module: General Study** (Type 2)
```cpp
gdcm::Attribute<0x0020, 0x000d> study_instance_uid;
study_instance_uid.Set(generate_uid().c_str());
ds.Insert(study_instance_uid.GetAsDataElement());

gdcm::Attribute<0x0008, 0x0020> study_date;
study_date.Set(get_current_date().c_str());
ds.Insert(study_date.GetAsDataElement());
```

**Module: RT Series** (Type 1 - Required)
```cpp
gdcm::Attribute<0x0008, 0x0060> modality;
modality.Set("RTDOSE");
ds.Insert(modality.GetAsDataElement());

gdcm::Attribute<0x0020, 0x000e> series_instance_uid;
series_instance_uid.Set(generate_uid().c_str());
ds.Insert(series_instance_uid.GetAsDataElement());
```

**Module: Frame of Reference** (Type 1)
```cpp
gdcm::Attribute<0x0020, 0x0052> frame_of_reference_uid;
frame_of_reference_uid.Set(generate_uid().c_str());
ds.Insert(frame_of_reference_uid.GetAsDataElement());
```

**Module: General Equipment** (Type 3 - Optional)
```cpp
gdcm::Attribute<0x0008, 0x0070> manufacturer;
manufacturer.Set("MOQUI");
ds.Insert(manufacturer.GetAsDataElement());

gdcm::Attribute<0x0018, 0x1020> software_versions;
software_versions.Set("1.0");
ds.Insert(software_versions.GetAsDataElement());
```

**Module: RT Dose** (Type 1)
```cpp
// Image dimensions
gdcm::Attribute<0x0028, 0x0010> rows;
rows.Set(ny);
ds.Insert(rows.GetAsDataElement());

gdcm::Attribute<0x0028, 0x0011> columns;
columns.Set(nx);
ds.Insert(columns.GetAsDataElement());

// Dose attributes
gdcm::Attribute<0x3004, 0x0001> dose_units;
dose_units.Set("GY");  // Gray
ds.Insert(dose_units.GetAsDataElement());

gdcm::Attribute<0x3004, 0x000e> dose_grid_scaling_attr;
dose_grid_scaling_attr.Set(dose_grid_scaling);
ds.Insert(dose_grid_scaling_attr.GetAsDataElement());

// ✅ FIXED: Spatial attributes now use correct separator strings
// All string formatting bugs have been fixed with proper escaped backslashes

// Pixel spacing - ✅ FIXED AT LINE 805
gdcm::Attribute<0x0028, 0x0030> pixel_spacing;
std::stringstream ps_ss;
ps_ss << std::fixed << std::setprecision(6) << dy << "\\" << dx;  // CORRECT
pixel_spacing.Set(ps_ss.str().c_str());
ds.Insert(pixel_spacing.GetAsDataElement());

// Grid frame offset vector - ✅ FIXED AT LINE 812
std::stringstream gf_ss;
for (int i = 0; i < nz; i++) {
    if (i > 0) gf_ss << "\\";  // CORRECT
    gf_ss << std::fixed << std::setprecision(6) << (grid.z_[i] - z0);
}
gdcm::Attribute<0x3004, 0x000c> grid_frame_offset_vector;
grid_frame_offset_vector.Set(gf_ss.str().c_str());
ds.Insert(grid_frame_offset_vector.GetAsDataElement());

// Image position (top-left corner of first voxel) - ✅ FIXED AT LINES 822-824
std::stringstream ip_ss;
ip_ss << std::fixed << std::setprecision(6)
      << (grid.x_[0] - dx/2) << "\\"  // CORRECT
      << (grid.y_[0] - dy/2) << "\\"  // CORRECT
      << (grid.z_[0] - dz/2);         // CORRECT
gdcm::Attribute<0x0020, 0x0032> image_position_patient;
image_position_patient.Set(ip_ss.str().c_str());
ds.Insert(image_position_patient.GetAsDataElement());

// Image orientation - ✅ FIXED AT LINE 831
gdcm::Attribute<0x0020, 0x0037> image_orientation_patient;
image_orientation_patient.Set("1\\0\\0\\0\\1\\0");  // CORRECT
ds.Insert(image_orientation_patient.GetAsDataElement());

// Pixel data
gdcm::Attribute<0x7fe0, 0x0010> pixel_data;
pixel_data.SetByteValue(reinterpret_cast<const char*>(pixel_data.data()),
                       length * sizeof(uint16_t));
ds.Insert(pixel_data.GetAsDataElement());
```

**✅ CRITICAL BUGS FIXED** in [mqi_io.hpp:805-831](moqui/base/mqi_io.hpp#L805-L831):
- Line 805: Fixed `dy << "\" << dx` to `dy << "\\" << dx`
- Line 812: Fixed `gf_ss << "\"` to `gf_ss << "\\"`
- Lines 822-824: Fixed `<< "\""` to `<< "\\"`
- Line 831: Fixed `"1\0\0\0\1\0"` to `"1\\0\\0\\0\\1\\0"`

All DICOM backslash separators are now correctly formatted!

### Step 5: Integration with save_reshaped_files() - ✅ COMPLETE

**Status**: DCM format handling is already integrated.

**Location**: [moqui/base/environments/mqi_tps_env.hpp:1611-1617](moqui/base/environments/mqi_tps_env.hpp#L1611-L1617)

**Actual integrated code**:
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
} else if (!this->output_format.compare("dcm")) {  // ✅ ALREADY IMPLEMENTED
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

### Step 6: Configuration Update - ✅ COMPLETE

**Status**: No changes needed - `OutputFormat` parameter already exists.

**Configuration file**: `tps_env/moqui_tps.in` line 42

**Usage**:
```
OutputFormat dcm
```

### Step 7: Helper Functions - ✅ COMPLETE

**Status**: All helper functions are implemented at [mqi_io.hpp:616-635](moqui/base/mqi_io.hpp#L616-L635)

**UID Generation** (uses GDCM's built-in UID generator):
```cpp
static std::string generate_uid() {
    gdcm::UIDGenerator uid_gen;
    return uid_gen.Generate();  // Uses GDCM's compliant UID generation
}
```

**Date/Time Functions** (implemented exactly as planned):
```cpp
static std::string get_current_date() {
    // Format: YYYYMMDD
    time_t now = time(0);
    struct tm tstruct = *localtime(&now);
    char buf[9];
    strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
    return buf;
}

static std::string get_current_time() {
    // Format: HHMMSS
    time_t now = time(0);
    struct tm tstruct = *localtime(&now);
    char buf[7];
    strftime(buf, sizeof(buf), "%H%M%S", &tstruct);
    return buf;
}
```

**Note**: Implementation uses GDCM's `UIDGenerator` instead of manual UID construction, which is more reliable and ensures DICOM compliance.

## Validation and Testing Strategy

**Purpose**: Validate the existing DCM export implementation and identify bugs.

### Test Case 1: 2cm Mode Output - CRITICAL BUG TESTING

**Configuration**: Use existing [2cmmode.md](2cmmode.md) setup
```
TwoCentimeterMode true
OutputFormat dcm
```

**Expected output**:
- File: `{beam_name}_1_{scorer_name}.dcm`
- Dimensions: 400×400×1 (Rows=400, Columns=400, NumberOfFrames=1)
- ⚠️ **BUG**: Pixel spacing will be malformed due to string bug at line 805
- ⚠️ **BUG**: Image position will be malformed due to bugs at lines 822-824
- ⚠️ **BUG**: Image orientation will contain null bytes instead of backslashes (line 831)
- Dose units: GY
- Data type: 16-bit unsigned integers

**Validation steps**:
1. **Generate DCM file** with current implementation
2. **DICOM conformance check**:
   ```bash
   dcmdump {beam_name}_1_{scorer_name}.dcm | grep -E "PixelSpacing|ImagePositionPatient|ImageOrientationPatient"
   ```
3. **Expected bug manifestations**:
   - PixelSpacing: Will show malformed separator (may show `"` instead of `\`)
   - ImagePositionPatient: Will show malformed separators
   - ImageOrientationPatient: Will show null characters or truncated values
4. **Viewer validation**:
   - Open in DICOM viewer (3D Slicer, MIM, RayStation)
   - ⚠️ May fail to load or show incorrect spatial positioning
   - Verify dose values are scaled correctly (if viewer can load file)
5. **DICOM validator**:
   - Run through official DICOM validator
   - Expected errors related to malformed multi-valued attributes

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

### Test Case 3: Compare with MHD Output (After Bug Fixes)

**Purpose**: Validate dose calculation correctness (independent of spatial metadata)

**Test**: Run same beam plan with both formats
```
OutputFormat mhd  # First run
OutputFormat dcm  # Second run (after fixing bugs)
```

**Validation**:
1. Load MHD in viewer (known working format)
2. Load DCM in same viewer (after bug fixes)
3. Subtract dose distributions
4. Verify difference is within numerical precision (< 0.001%)
5. Verify spatial alignment is identical

### Test Case 4: DICOM Conformance Testing

**Tools required**:
- `dcmdump` (part of DCMTK toolkit): For inspecting DICOM tags generated by GDCM
- `dciodvfy` (DICOM validator): For IOD conformance checking
- 3D Slicer or similar: For visual validation

**Validation commands**:
```bash
# Dump DICOM structure
dcmdump output.dcm > output_dump.txt

# Check critical attributes
dcmdump output.dcm | grep -E "(0028,0030|0020,0032|0020,0037|3004,000c)"

# Validate IOD conformance
dciodvfy output.dcm

# Check for warnings/errors in spatial attributes
```

**Expected findings (before bug fix)**:
- Errors in multi-valued string attributes
- Warnings about separator characters
- Possible rejection by strict DICOM validators

## Implementation Status Checklist

### Completed Implementation ✅
- [x] **CMake**: GDCM dependency integrated ([CMakeLists.txt:37-68](tps_env/CMakeLists.txt#L37-L68))
- [x] **mqi_io.hpp**: `save_to_dcm()` declaration at line 31 ([link](moqui/base/mqi_io.hpp#L31))
- [x] **mqi_io.hpp**: `save_to_dcm()` implementation at lines 637-856 ([link](moqui/base/mqi_io.hpp#L637-L856))
- [x] **mqi_io.hpp**: Helper functions at lines 616-635 ([link](moqui/base/mqi_io.hpp#L616-L635))
- [x] **mqi_tps_env.hpp**: DCM format handling integrated at lines 1611-1617 ([link](moqui/base/environments/mqi_tps_env.hpp#L1611-L1617))
- [x] **Build**: Compiles successfully with GDCM

### Critical Bug Fixes Completed ✅
- [x] **mqi_io.hpp:805**: Fixed pixel spacing separator - Changed `"\"` to `"\\"`
- [x] **mqi_io.hpp:812**: Fixed grid frame offset separator - Changed `"\";` to `"\\";`
- [x] **mqi_io.hpp:822-824**: Fixed image position separators - Changed `"\""` to `"\\"`
- [x] **mqi_io.hpp:831**: Fixed image orientation - Changed `"1\0\0\0\1\0"` to `"1\\0\\0\\0\\1\\0"`

### Testing and Validation Required 🧪
- [ ] **Test**: Generate DCM file in 2cm mode and inspect with dcmdump
- [ ] **Test**: Validate bug manifestations with DICOM viewer
- [ ] **Test**: Regular phantom with DCM output
- [ ] **Test**: Compare DCM vs MHD dose distributions (after bug fixes)
- [ ] **Test**: DICOM conformance validation with dciodvfy
- [ ] **Test**: Load in clinical treatment planning system

### Documentation
- [x] **plan.md**: Updated to reflect actual implementation status
- [ ] **bug_report.md**: Create dedicated bug fix documentation
- [ ] **user_guide**: Update with DCM format usage and known issues

## File Summary

| File | Lines | Implementation Status | Bug Fixes Required |
|------|-------|----------------------|-------------------|
| `tps_env/CMakeLists.txt` | 37-68 | ✅ Complete (GDCM integrated) | None |
| `moqui/base/mqi_io.hpp` | 31-38 | ✅ Complete (declaration) | None |
| `moqui/base/mqi_io.hpp` | 637-856 | ✅ Complete (220 lines) | **All bugs fixed** |
| `moqui/base/mqi_io.hpp` | 616-635 | ✅ Complete (helper functions) | None |
| `moqui/base/environments/mqi_tps_env.hpp` | 1611-1617 | ✅ Complete (DCM handling) | None |
| `tps_env/moqui_tps.in` | Line 42 | ✅ Complete (OutputFormat param) | None |

**Critical bugs fixed**:
1. Line 805: Pixel spacing separator ✅
2. Line 812: Grid frame offset separator ✅
3. Lines 822-824: Image position separators ✅
4. Line 831: Image orientation string ✅

## Related Documentation

- **Inconsistencies Report**: See [plan_inconsistencies_report.md](plan_inconsistencies_report.md) for detailed analysis
- **2cm Mode Details**: See [2cmmode.md](2cmmode.md) for phantom geometry and scorer configuration
- **DICOM Standard**: DICOM PS3.3 Section C.8.8 (RT Dose Module)
- **GDCM Documentation**: https://gdcm.sourceforge.net/html/index.html
- **GDCM API Reference**: https://gdcm.sourceforge.net/html/annotated.html
- **Current Implementation**: [moqui/base/mqi_io.hpp:637-856](moqui/base/mqi_io.hpp#L637-L856)

## Summary

This document has been updated to reflect the **actual completed state** of the DICOM RT Dose export implementation:

1. ✅ **Implementation is complete** - Full DCM export functionality exists
2. ✅ **GDCM library is used** - Not DCMTK as originally planned
3. ✅ **Critical bugs fixed** - All 4 string formatting bugs have been resolved
4. 🧪 **Testing required** - Validation and verification needed

**Next steps**:
1. ✅ Fixed the 4 critical string formatting bugs in [mqi_io.hpp:805-831](moqui/base/mqi_io.hpp#L805-L831)
2. Validate DCM output with DICOM conformance tools
3. Test with clinical DICOM viewers and treatment planning systems

## Bug Fix Documentation

### Critical String Formatting Bugs (FIXED)

Four critical string formatting bugs were identified and fixed in the DICOM RT Dose export implementation. These bugs caused incorrect separator characters to be written to DICOM files, which would corrupt spatial attributes and prevent proper loading in DICOM viewers.

#### Bug Details

**Location**: [`moqui/base/mqi_io.hpp:805-831`](moqui/base/mqi_io.hpp#L805-L831)

**Root Cause**: Incorrect escape sequences for backslash characters in DICOM multi-valued strings. DICOM requires backslash (`\`) as a separator between values, but the code was using incorrect escape sequences.

#### Fixed Issues

1. **Line 805 - Pixel Spacing Attribute**
   - **Buggy Code**: `ps_ss << std::fixed << std::setprecision(6) << dy << "\" << dx;`
   - **Fixed Code**: `ps_ss << std::fixed << std::setprecision(6) << dy << "\\" << dx;`
   - **Impact**: Pixel spacing values were malformed, causing incorrect voxel size interpretation

2. **Line 812 - Grid Frame Offset Vector**
   - **Buggy Code**: `if (i > 0) gf_ss << "\";`
   - **Fixed Code**: `if (i > 0) gf_ss << "\\";`
   - **Impact**: Z-axis position values were malformed, causing incorrect slice positioning

3. **Lines 822-824 - Image Position Patient**
   - **Buggy Code**:
     ```
     << (grid.x_[0] - dx/2) << "\"
     << (grid.y_[0] - dy/2) << "\"
     << (grid.z_[0] - dz/2);
     ```
   - **Fixed Code**:
     ```
     << (grid.x_[0] - dx/2) << "\\"
     << (grid.y_[0] - dy/2) << "\\"
     << (grid.z_[0] - dz/2);
     ```
   - **Impact**: 3D position of the image was malformed, causing incorrect spatial registration

4. **Line 831 - Image Orientation Patient**
   - **Buggy Code**: `image_orientation_patient.Set("1\0\0\0\1\0");`
   - **Fixed Code**: `image_orientation_patient.Set("1\\0\\0\\0\\1\\0");`
   - **Impact**: Image orientation vectors contained null bytes instead of backslashes, causing incorrect orientation

#### Validation

After fixes, DICOM files should have properly formatted spatial attributes:
- PixelSpacing: "1.000000\1.000000" (example)
- ImagePositionPatient: "-200.000000\-200.000000\-1.000000" (example)
- ImageOrientationPatient: "1\0\0\0\1\0" (correct format)

#### Testing Recommendation

To verify the fixes:
1. Generate a DICOM RT Dose file with `OutputFormat dcm`
2. Use `dcmdump` to inspect the spatial attributes:
   ```bash
   dcmdump output.dcm | grep -E "(0028,0030|0020,0032|0020,0037)"
   ```
3. Verify that backslashes appear as separators (not quotes or null bytes)
4. Load the DICOM file in a viewer to confirm proper spatial positioning

## Troubleshooting Guide

### Common Issues and Solutions

#### DICOM File Won't Open in Viewer

**Symptoms**:
- DICOM viewer reports "invalid file format" or "corrupted file"
- Spatial attributes appear malformed
- Images fail to load or display incorrectly

**Possible Causes**:
1. **String formatting bugs** (FIXED in current implementation)
   - Check for malformed backslash separators in spatial attributes
   - Use `dcmdump` to inspect attributes: `dcmdump file.dcm | grep -E "(0028,0030|0020,0032|0020,0037)"`

2. **Missing required DICOM attributes**
   - Verify all required RT Dose attributes are present
   - Check with DICOM validator: `dciodvfy file.dcm`

3. **Incorrect transfer syntax**
   - Ensure Explicit VR Little Endian is used
   - Check with: `dcmdump file.dcm | grep TransferSyntaxUID`

#### Incorrect Dose Values

**Symptoms**:
- Dose values appear as 0 or extremely large/small
- Dose doesn't match expected distribution

**Solutions**:
1. **Check dose grid scaling**
   - Verify DoseGridScaling attribute is reasonable (typically 1e-6 to 1e-3)
   - Inspect with: `dcmdump file.dcm | grep DoseGridScaling`

2. **Verify scaling factor**
   - Check `particles_per_history` value in configuration
   - Ensure dose is being scaled correctly before DICOM conversion

3. **Compare with other formats**
   - Generate same plan with MHD format
   - Compare dose values using external tools

#### Spatial Positioning Issues

**Symptoms**:
- Dose doesn't overlay correctly on CT anatomy
- Incorrect coordinate system or orientation
- Voxel spacing appears wrong

**Solutions**:
1. **Verify coordinate system**
   - Check ImagePositionPatient and ImageOrientationPatient
   - Ensure values match expected phantom geometry

2. **Validate voxel spacing**
   - Check PixelSpacing attribute matches expected voxel size
   - For 2cm mode: should be approximately "1.0\1.0" (mm)

3. **Check grid dimensions**
   - Verify Rows, Columns, and NumberOfFrames match grid size
   - For 2cm mode: should be 400×400×1

#### Build/Compilation Issues

**Symptoms**:
- GDCM headers not found during compilation
- Linker errors for GDCM libraries

**Solutions**:
1. **Verify GDCM installation**
   - Ensure GDCM is properly installed and in PATH
   - Check CMakeLists.txt includes correct GDCM modules

2. **Check include paths**
   - Verify GDCM include directories are accessible
   - Update CMake include_directories if needed

3. **Validate library linking**
   - Ensure all required GDCM libraries are linked
   - Check for missing gdcmRT module specifically

### DICOM Conformance Validation

#### Required Tools
- `dcmdump` (DCMTK): For inspecting DICOM tags
- `dciodvfy` (DCMTK): For IOD conformance checking
- DICOM viewer (3D Slicer, MIM, etc.): For visual validation

#### Validation Commands
```bash
# Basic structure validation
dcmdump output.dcm > output_dump.txt

# Check critical attributes
dcmdump output.dcm | grep -E "(0028,0030|0020,0032|0020,0037|3004,000c|3004,000e)"

# Full IOD conformance check
dciodvfy output.dcm

# Check for warnings/errors
dciodvfy -v output.dcm
```

#### Expected Validation Results
After bug fixes, the DICOM file should:
- Pass dciodvfy validation without errors
- Show properly formatted multi-valued attributes
- Load correctly in standard DICOM viewers
- Display dose with correct spatial positioning

### Performance Considerations

#### Large Dose Arrays
- For large 3D dose arrays, consider memory usage during conversion
- Monitor memory usage when converting from double to uint16_t
- Consider chunked processing for very large datasets

#### File Size Optimization
- DICOM RT Dose uses 16-bit unsigned integers for pixel data
- This provides good balance between precision and file size
- DoseGridScaling attribute ensures proper dose representation

## Known Limitations

### Current Implementation Constraints

#### DICOM Compliance
- **Limited RT Module Support**: Only RT Dose Storage is implemented
  - RT Plan, RT Structure Set, and RT Image not supported
  - No support for dose references or treatment plan relationships
- **Simplified Patient Information**: Uses placeholder values
  - Patient Name: "PHANTOM"
  - Patient ID: "QA_PHANTOM"
  - No real patient demographics integration

#### Spatial Accuracy
- **Coordinate System Assumptions**:
  - Assumes standard HFS (Head First Supine) orientation
  - No support for patient reorientation (LPO, RPO, etc.)
  - Fixed image orientation: "1\0\0\0\1\0"
- **Grid Limitations**:
  - Only supports regular, rectangular grids
  - No support for non-uniform voxel spacing
  - Limited to single-frame or simple multi-frame dose distributions

#### Data Representation
- **Dose Precision**:
  - 16-bit unsigned integer representation limits precision
  - Very small dose values may be lost due to scaling
  - DoseGridScaling may not preserve very low-dose regions
- **Dynamic Range**:
  - Single DoseGridScaling factor for entire dataset
  - May not optimize for both high and low dose regions simultaneously

#### Integration Limitations
- **No Referenced Objects**:
  - No links to referenced RT Plan or Structure Set
  - No support for dose ROI mapping
  - No beam sequence information
- **Missing Metadata**:
  - No treatment machine information
  - No beam energy or modality details
  - No treatment date or fraction information

#### Performance Considerations
- **Memory Usage**:
  - Creates full dense array before DICOM conversion
  - May be inefficient for very sparse dose distributions
  - No streaming or chunked processing for large datasets
- **File Size**:
  - No built-in compression for pixel data
  - All frames stored in single file (no multi-file optimization)

### Potential Enhancements

#### DICOM Extensions
1. **RT Plan Integration**
   - Add support for referenced RT Plan
   - Include beam sequence and fraction information
   - Support for multiple treatment phases

2. **Structure Set Support**
   - Implement RT Structure Set export
   - Enable dose-volume histogram (DVH) calculation
   - Support for ROI-based dose analysis

3. **Enhanced Metadata**
   - Treatment machine specifications
   - Beam quality and energy details
   - Treatment session information

#### Technical Improvements
1. **Compression Support**
   - Add DICOM compression options
   - Support for lossless and lossy compression
   - Optimized storage for large datasets

2. **Coordinate System Flexibility**
   - Support for patient reorientation
   - Flexible image orientation matrices
   - Non-standard coordinate systems

3. **Performance Optimization**
   - Streaming processing for large arrays
   - Memory-efficient sparse-to-dense conversion
   - Parallel processing for multi-core systems

### Clinical Considerations

#### QA Workflow Integration
- **No Automatic QA Checks**: No built-in validation of dose distribution
- **Missing Comparison Tools**: No support for gamma analysis or dose comparison
- **Limited Reporting**: No automatic report generation

#### Regulatory Compliance
- **Audit Trail**: No logging of dose calculation parameters
- **Version Control**: Limited tracking of calculation versions
- **Access Control**: No built-in user authentication or authorization

### Platform Dependencies

#### Build Requirements
- **GDCM Library**: Requires specific version compatibility
- **C++ Compiler**: Needs C++11 or later for some features
- **System Dependencies**:
  - `sys/mman.h` for memory mapping (Unix/Linux specific)
  - May need Windows equivalents for cross-platform support

#### Runtime Dependencies
- **DICOM Viewers**: Requires compatible DICOM viewer for validation
- **Validation Tools**: DCMTK tools needed for conformance checking
- **Treatment Planning Systems**: May need specific import capabilities