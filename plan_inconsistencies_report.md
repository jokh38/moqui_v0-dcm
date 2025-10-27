# Plan.md Inconsistencies Report

## Executive Summary

The [plan.md](plan.md) file contains a comprehensive implementation plan for DICOM RT Dose export, **but the implementation is already complete**. This document identifies 10 specific inconsistencies between the plan and the actual codebase that need correction.

**Critical Finding**: The plan narrative treats DCM export as "to be implemented" when it's actually **already fully implemented and integrated** (lines 637-856 in [mqi_io.hpp](moqui/base/mqi_io.hpp:637-856)).

---

## Detailed Inconsistencies

### 1. Library Mismatch: DCMTK vs GDCM

**Plan states** (multiple locations):
- Line 205: "Add DCMTK Dependency to Build System"
- Lines 212-214: References DCMTK package finding
- Line 218: Lists "dcmdata" and "dcmrt" DCMTK modules

**Actual implementation**:
- [mqi_io.hpp:15-20](moqui/base/mqi_io.hpp:15-20): Uses GDCM headers
  ```cpp
  #include "gdcmAttribute.h"
  #include "gdcmFile.h"
  #include "gdcmWriter.h"
  #include "gdcmUIDGenerator.h"
  #include "gdcmMediaStorage.h"
  ```
- [CMakeLists.txt:37-38](tps_env/CMakeLists.txt:37-38): `find_package(GDCM REQUIRED)`
- [CMakeLists.txt:55-68](tps_env/CMakeLists.txt:55-68): Links GDCM libraries (gdcmCommon, gdcmDICT, gdcmIOD, etc.)

**Fix required**: Replace all references to DCMTK with GDCM throughout the plan.

---

### 2. Function Signature Type Mismatch

**Plan shows** (lines 227-233):
```cpp
template<typename R>
void save_to_dcm(const mqi::node<R>*  child,
                 const double*        src,
                 // ...
```

**Actual implementation** ([mqi_io.hpp:31-38](moqui/base/mqi_io.hpp:31-38)):
```cpp
template<typename R>
void
save_to_dcm(const mqi::node_t<R>* children,  // ← Note: node_t, not node
            const double*         src,
            // ...
```

**Fix required**: Change `mqi::node<R>*` to `mqi::node_t<R>*` throughout the plan.

---

### 3. Grid Geometry Access Pattern

**Plan shows** (line 262):
```cpp
const mqi::grid3d<R>& grid = child->geo_;
```

**Actual implementation** ([mqi_io.hpp:647](moqui/base/mqi_io.hpp:647)):
```cpp
const mqi::grid3d<R>& grid = children->geo_[0];  // ← Array access [0]
```

**Fix required**: Update to show `children->geo_[0]` to reflect that geo_ is an array.

---

### 4. String Formatting Errors (Critical Bug Documentation)

**Plan shows incorrect escape sequences** at multiple locations:

| Line | Plan Code | Should Be | Issue |
|------|-----------|-----------|-------|
| 805 | `ps_ss << dy << "\" << dx;` | `ps_ss << dy << "\\" << dx;` | Missing backslash escape |
| 812 | `if (i > 0) gf_ss << "\";` | `if (i > 0) gf_ss << "\\";` | Missing backslash escape |
| 822-824 | `<< (grid.x_[0] - dx/2) << "\"`<br>`<< (grid.y_[0] - dy/2) << "\"`<br>`<< (grid.z_[0] - dz/2);` | `<< (grid.x_[0] - dx/2) << "\\"`<br>`<< (grid.y_[0] - dy/2) << "\\"`<br>`<< (grid.z_[0] - dz/2);` | Missing backslash escapes |
| 831 | `image_orientation_patient.Set("1\0\0\0\1\0");` | `image_orientation_patient.Set("1\\0\\0\\0\\1\\0");` | Incorrect null characters |

**Actual implementation** ([mqi_io.hpp:805-831](moqui/base/mqi_io.hpp:805-831)) has the **exact same bugs**:
```cpp
// Line 805
ps_ss << std::fixed << std::setprecision(6) << dy << "\" << dx;

// Line 812
if (i > 0) gf_ss << "\";

// Lines 822-824
ip_ss << std::fixed << std::setprecision(6)
      << (grid.x_[0] - dx/2) << "\"
      << (grid.y_[0] - dy/2) << "\"
      << (grid.z_[0] - dz/2);

// Line 831
image_orientation_patient.Set("1\0\0\0\1\0");
```

**Impact**: These bugs will cause **DICOM file corruption**:
- Pixel spacing, grid frame offsets, and image position will have incorrect separators
- DICOM backslash separators (`\`) won't be properly written
- Image orientation will contain null characters instead of backslashes

**Fix required**:
1. Correct the string formatting in plan.md
2. **Create separate bug report for actual implementation bugs**
3. Fix the actual code in [mqi_io.hpp:805-831](moqui/base/mqi_io.hpp:805-831)

---

### 5. Implementation Status Misrepresentation

**Plan narrative** (Overview, line 3-5):
> "This document describes the **implementation plan** for **adding** DICOM RT Dose export capability to MOQUI."

**Plan narrative** (Line 203):
> "## DICOM RT Dose Export Implementation Plan"

**Plan checklist** (Lines 518-527):
```
- [ ] **mqi_io.hpp**: Add `save_to_dcm()` function declaration (after line 98)
- [ ] **mqi_io.hpp**: Implement `save_to_dcm()` function (after line 597)
- [ ] **mqi_io.hpp**: Add helper functions (UID generation, date/time)
```

**Actual status**:
- DCM export is **fully implemented** ([mqi_io.hpp:637-856](moqui/base/mqi_io.hpp:637-856))
- Helper functions exist ([mqi_io.hpp:616-635](moqui/base/mqi_io.hpp:616-635))
- Integration is complete ([mqi_tps_env.hpp:1611-1617](moqui/base/environments/mqi_tps_env.hpp:1611-1617))

**Fix required**:
- Change title from "Implementation Plan" to "Implementation Validation Guide"
- Reframe narrative from "adding" to "validating and improving"
- Update checklist from "implement" to "verify" or "test"

---

### 6. Build System Status

**Plan states** (Lines 205-215):
> "### Step 1: Add DCMTK Dependency to Build System
>
> **Files to modify**:
> - `CMakeLists.txt` - Add DCMTK package finding and linking"

**Actual status** ([CMakeLists.txt:37-68](tps_env/CMakeLists.txt:37-68)):
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
    gdcmRT
)
```

**Fix required**: Update Step 1 to reflect that GDCM is already integrated, change from "Add" to "Verify".

---

### 7. Configuration Integration

**Plan states** (Lines 382-416):
> "### Step 5: Integration with save_reshaped_files()
>
> **Location**: `moqui/base/environments/mqi_tps_env.hpp` lines 1597-1617
>
> **Add DCM format handling**: [shows else-if block to add]"

**Actual status** ([mqi_tps_env.hpp:1611-1617](moqui/base/environments/mqi_tps_env.hpp:1611-1617)):
```cpp
} else if (!this->output_format.compare("dcm")) {
    mqi::io::save_to_dcm<R>(this->world->children[c_ind],
                            reshaped_data,
                            this->particles_per_history,
                            this->output_path,
                            filename,
                            vol_size);
```

**Fix required**: Change Step 5 from "Integration with" to "Verification of integration in" since it's already done.

---

### 8. Function Declaration Location

**Plan states** (Lines 236-246):
> "**Function declaration** (add after line 98):
> ```cpp
> /// Save dose data as DICOM RT Dose file
> template<typename R>
> void save_to_dcm(const mqi::node<R>*  child,
> // ...
> ```"

**Actual location** ([mqi_io.hpp:30-38](moqui/base/mqi_io.hpp:30-38)):
- Declaration is at **line 31**, not after line 98
- It's in the global `mqi::` namespace, before the `mqi::io::` namespace begins

**Fix required**: Update line reference from "after line 98" to "at line 31".

---

### 9. Helper Functions Status

**Plan states** (Lines 427-459):
> "### Step 7: Helper Functions to Implement
>
> **UID Generation** (DICOM requires unique identifiers):
> ```cpp
> std::string generate_uid() {
>     // Format: root.timestamp.random
>     // ...
> ```"

**Actual status** ([mqi_io.hpp:616-635](moqui/base/mqi_io.hpp:616-635)):
```cpp
// Helper functions for DICOM RT Dose export
static std::string generate_uid() {
    gdcm::UIDGenerator uid_gen;
    return uid_gen.Generate();
}

static std::string get_current_date() {
    time_t now = time(0);
    struct tm tstruct = *localtime(&now);
    char buf[9];
    strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
    return buf;
}

static std::string get_current_time() {
    time_t now = time(0);
    struct tm tstruct = *localtime(&now);
    char buf[7];
    strftime(buf, sizeof(buf), "%H%M%S", &tstruct);
    return buf;
}
```

**Fix required**:
- Change "Helper Functions to Implement" to "Helper Functions (Already Implemented)"
- Update UID generation implementation to show GDCM's UIDGenerator usage instead of custom implementation

---

### 10. Testing Focus

**Plan states** (Lines 461-515):
> "## Testing Strategy
>
> ### Test Case 1: 2cm Mode Output
>
> **Expected output**:
> - File: `{beam_name}_1_{scorer_name}.dcm`
> - Dimensions: 400×400×1"

**Appropriate focus** (given implementation is complete):
- Should emphasize **validation** of existing implementation
- Should focus on **bug detection** (especially string formatting issues)
- Should include **DICOM conformance testing**
- Should add **comparison with reference implementations**

**Fix required**:
- Reframe from "Testing Strategy for New Feature" to "Validation Strategy for Existing Implementation"
- Add test cases specifically for string separator bugs
- Add DICOM validator tool recommendations (dcmdump, DICOM validator)

---

## Implementation Checklist Status

**Plan checklist** (Lines 518-527) shows all items unchecked:
```
- [ ] **CMake**: Add DCMTK dependency
- [ ] **mqi_io.hpp**: Add `save_to_dcm()` function declaration (after line 98)
- [ ] **mqi_io.hpp**: Implement `save_to_dcm()` function (after line 597)
- [ ] **mqi_io.hpp**: Add helper functions (UID generation, date/time)
- [ ] **mqi_tps_env.hpp**: Add DCM format handling in `save_reshaped_files()`
```

**Actual status**: All items are complete ✅

**Proposed revised checklist**:
```
- [x] **CMake**: GDCM dependency integrated (lines 37-68)
- [x] **mqi_io.hpp**: `save_to_dcm()` declaration at line 31
- [x] **mqi_io.hpp**: `save_to_dcm()` implementation at lines 637-856
- [x] **mqi_io.hpp**: Helper functions at lines 616-635
- [x] **mqi_tps_env.hpp**: DCM format handling at lines 1611-1617
- [ ] **BUG FIX REQUIRED**: String formatting errors at lines 805, 812, 822-824, 831
- [ ] **TEST**: Validate 2cm mode DCM output
- [ ] **TEST**: Validate regular phantom DCM output
- [ ] **TEST**: Compare DCM vs MHD dose distributions
- [ ] **TEST**: DICOM conformance validation with dcmdump
```

---

## Recommended Actions

### Immediate (High Priority)

1. **Fix Critical Bugs in Implementation**
   - File: [mqi_io.hpp:805-831](moqui/base/mqi_io.hpp:805-831)
   - Issue: String formatting with incorrect escape sequences
   - Impact: DICOM file corruption
   - Action: Replace `"\"` with `"\\"` and fix image orientation string

2. **Update Plan Narrative**
   - Change from "Implementation Plan" to "Implementation Validation and Bug Fix Guide"
   - Reframe all "to be implemented" language to "to be validated/tested"

### Secondary (Documentation)

3. **Correct Technical Details**
   - Replace DCMTK → GDCM throughout
   - Fix function signatures (node<R> → node_t<R>)
   - Fix grid access (geo_ → geo_[0])
   - Update line number references

4. **Update Testing Section**
   - Add validation focus
   - Add DICOM conformance testing steps
   - Add bug detection test cases

### Tertiary (Enhancement)

5. **Add Missing Documentation**
   - Document GDCM API usage patterns
   - Add troubleshooting section
   - Add known limitations section

---

## File Modification Summary

| File | Current State | Required Action |
|------|---------------|-----------------|
| [plan.md](plan.md) | Describes unimplemented feature | Update narrative: implementation → validation |
| [mqi_io.hpp:805-831](moqui/base/mqi_io.hpp:805-831) | Contains string formatting bugs | **FIX BUGS**: Correct escape sequences |
| [mqi_io.hpp:637-856](moqui/base/mqi_io.hpp:637-856) | Complete implementation | Document and validate |
| [CMakeLists.txt:37-68](tps_env/CMakeLists.txt:37-68) | GDCM already integrated | Update plan to reflect |
| [mqi_tps_env.hpp:1611-1617](moqui/base/environments/mqi_tps_env.hpp:1611-1617) | DCM handling integrated | Update plan to reflect |

---

## Conclusion

The plan.md document provides excellent technical documentation but is **fundamentally misaligned** with the actual project state. The DCM export feature is **not planned—it's implemented**. The document needs a comprehensive revision to:

1. **Fix critical bugs** identified in the string formatting
2. **Shift focus** from implementation to validation
3. **Correct technical details** (DCMTK→GDCM, type names, line numbers)
4. **Update status** of all completed items

**Recommended next step**: Create a separate "DCM Export Bug Fix Plan" to address the string formatting issues before they cause production problems.
