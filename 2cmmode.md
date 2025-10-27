# 2cm Mode Analysis - MOQUI TPS

## Overview

The 2cm mode (`TwoCentimeterMode true` in `moqui_tps.in`) is designed for patient-specific QA dose measurements using a simplified water phantom geometry. This document describes the actual implementation based on code analysis.

## Configuration

**Location**: `moqui_tps.in` line 17
```
TwoCentimeterMode true
```

**Implementation**: `moqui/base/environments/mqi_tps_env.hpp` lines 787-986

## Beam Configuration

### Beam Direction
When `usingPhantomGeo` is enabled (lines 823-833, 1088-1099):
- All beam angles are set to **0 degrees**:
  ```cpp
  p_coord.angles[0] = 0.f;  // Gantry angle
  p_coord.angles[1] = 0.f;  // Couch angle
  p_coord.angles[2] = 0.f;
  p_coord.angles[3] = 0.f;  // IEC2DICOM angle
  ```
- Coordinate translations are zeroed
- **Beam travels in NEGATIVE Z direction** (from snout toward isocenter)

### Air Box (lines 843-862)
Creates air gap between snout and phantom:
```cpp
X: [-200, 200] mm  (1 voxel)
Y: [-200, 200] mm  (1 voxel)
Z: [20, snoutPos] mm  (1 voxel, ~300-400mm depending on snout position)
```
Filled with air density.

## Phantom Geometry

The 2cm mode creates **three sequential water phantom slabs** arranged along the Z-axis:

### Overall Layout
```
Beam Direction (−Z) →

[Air Box]  →  [Front Phantom]  →  [Center Phantom]  →  [Back Phantom]
   ↓              ↓                     ↓                     ↓
+400mm         +20 to +1            +1 to −1              −1 to −380
  (air)       (19mm water)        (2mm water)           (379mm water)
               NO SCORER           ★ SCORER ★             NO SCORER
```

### Front Phantom (lines 930-947)
**Purpose**: Dose buildup region before measurement plane

```cpp
grid3d(-200, 200, 401,    // X: 400 voxels, 1mm each
       -200, 200, 401,    // Y: 400 voxels, 1mm each
       1, 20, 2)          // Z: 1 voxel, 19mm thick
```
- Dimensions: **400×400×1 voxels**
- Physical size: **400×400×19 mm³**
- Material: Water (ρ = 1.0 g/cm³)
- **No scorer attached**

### Center Phantom (lines 949-966) - THE SCORER
**Purpose**: Measurement plane for 2D dose distribution

```cpp
grid3d(-200, 200, 401,    // X: 400 voxels, 1mm each
       -200, 200, 401,    // Y: 400 voxels, 1mm each
       -1, 1, 2)          // Z: 1 voxel, 2mm thick
```
- Dimensions: **400×400×1 voxels**
- Physical size: **400×400×2 mm³**
- Voxel size: **1mm × 1mm × 2mm**
- Material: Water (ρ = 1.0 g/cm³)
- **Scorer IS attached** (lines 1037-1050)
- Position: Centered at Z ≈ 0 (±1mm)

### Back Phantom (lines 969-986)
**Purpose**: Backscatter region after measurement plane

```cpp
grid3d(-200, 200, 401,    // X: 400 voxels, 1mm each
       -200, 200, 401,    // Y: 400 voxels, 1mm each
       -380, -1, 2)       // Z: 1 voxel, 379mm thick
```
- Dimensions: **400×400×1 voxels**
- Physical size: **400×400×379 mm³**
- Material: Water (ρ = 1.0 g/cm³)
- **No scorer attached**

## Scorer Configuration

### Scorer Attachment (lines 1022-1050)
Only the **center phantom** node receives a scorer:
```cpp
phantom->n_scorers = 1;
phantom->scorers[0] = new mqi::scorer<R>(
    this->scorer_string.c_str(),
    this->dcm_.dim_.x * this->dcm_.dim_.y * this->dcm_.dim_.z,  // 400×1×400 = 160,000
    fp0
);
```

### Scorer Dimensions
- **Grid**: 400×400×1 voxels
- **Total voxels**: 160,000
- **Output**: 2D dose map in XY plane
- **Plane orientation**: Perpendicular to beam direction (Z-axis)

### What Gets Recorded
The scorer records dose deposition in a **400×400 2D grid** representing:
- **X dimension**: 400 pixels, 1mm spacing (−200 to +200 mm)
- **Y dimension**: 400 pixels, 1mm spacing (−200 to +200 mm)
- **Z dimension**: 1 pixel, 2mm thick (−1 to +1 mm)

This creates a 2D dose distribution in the XY plane at the isocenter (Z ≈ 0).

## Output (lines 1583-1621)

The code iterates through all world children and saves only nodes with scorers:
```cpp
for (int c_ind = 0; c_ind < this->world->n_children; c_ind++) {
    for (int s_ind = 0; s_ind < this->world->children[c_ind]->n_scorers; s_ind++) {
        filename = beam_name + "_" + std::to_string(c_ind) + "_" + scorer_name;
        // Save dose data
    }
}
```

Since only the center phantom has a scorer, the output will be:
- **One file per beam**
- **Dimensions**: 400×400×1
- **Format**: Binary (.raw), MetaImage (.mhd/.mha), or NumPy (.npz)
- **Content**: 2D dose distribution in water at 2cm depth

## Physical Interpretation

### Why Three Phantoms?

1. **Front Phantom (19mm)**: Provides dose buildup region
   - Proton beams build up dose as they enter tissue
   - This region simulates initial buildup before measurement

2. **Center Phantom (2mm)**: The actual measurement plane
   - Thin slice minimizes dose averaging in beam direction
   - Located at isocenter for accurate geometry
   - Records the 2D dose distribution

3. **Back Phantom (379mm)**: Provides backscatter
   - Simulates tissue beyond measurement plane
   - Contributes realistic backscatter to dose at measurement depth
   - Total depth: 19 + 2 + 379 = 400mm = 40cm (sufficient for proton stopping)

### Coordinate System Summary

| Axis | Range | Resolution | Purpose |
|------|-------|------------|---------|
| **X** | −200 to +200 mm | 1mm (400 voxels) | Lateral dimension 1 |
| **Y** | −200 to +200 mm | 1mm (400 voxels) | Lateral dimension 2 |
| **Z** | −1 to +1 mm | 2mm (1 voxel) | Beam direction (measurement thickness) |

### Verification

✅ **2D dose recording confirmed**: The scorer is a 400×400×1 grid
✅ **Orthogonal to beam**: Scorer plane (XY) ⊥ beam direction (Z)
✅ **Correct positioning**: Centered at isocenter (Z ≈ 0)
✅ **Appropriate resolution**: 1mm lateral, 2mm longitudinal

## Key Implementation Details

### Variable Confusion and Coordinate System Transformation

The variable naming in lines 787-797 can be confusing because it involves **two different coordinate systems**:

```cpp
this->phantomDimX = 400;  // DICOM X dimension → World X-axis (400 voxels)
this->phantomDimY = 1;    // DICOM Y dimension → World Z-axis (1 voxel, 2mm thick)
this->phantomDimZ = 400;  // DICOM Z dimension → World Y-axis (400 voxels)
```

#### Coordinate System Mapping

There are two coordinate systems at play:

1. **DICOM Coordinate System** (used for `dcm_.dim_` and variable naming)
2. **World/Physics Coordinate System** (used for actual geometry construction)

The transformation is: **DICOM (X,Y,Z) → World (X,Z,Y)**

#### Why This Mapping Exists

1. **DICOM Convention**: In medical imaging, Z typically represents the slice direction (axial direction)
2. **Physics Convention**: In particle transport, Z typically represents the beam direction
3. **2cm Mode Requirements**: Creates a thin measurement plane perpendicular to the beam, so the "thin" dimension (1 voxel) needs to align with the beam direction

#### Evidence in Code

- **Line 458**: `dcm.dim_ = { this->phantomDimX, this->phantomDimY, this->phantomDimZ }` (DICOM coordinates)
- **Lines 950-959**: Actual grid construction uses transformed coordinates where scorer plane is XY in world coordinates
- **Line 1039**: Scorer capacity uses DICOM dimensions (400 * 1 * 400 = 160,000)

This explains why `phantomDimY` has "Y" in its name (following DICOM conventions) but the actual scorer plane is correctly defined in the XY direction in world coordinates through this coordinate system transformation.

The **actual phantom geometry** is defined separately in the grid3d constructors (lines 931-979) using the transformed world coordinate system.

### Node Hierarchy
```
world
├── beamline_objects[0..n]  (range shifters, apertures, etc.)
├── airBox                  (no scorer)
├── frontPhantom            (no scorer)
├── phantom                 (★ HAS SCORER ★)
└── backPhantom             (no scorer)
```

Only `phantom` (center slice) contributes to dose output.

## Usage Notes

When analyzing 2cm mode output:
- The output is a **2D array** (400×400) despite being stored as 3D (400×400×1)
- The dose represents measurement at **isocenter plane** (Z = 0)
- The measurement includes realistic **buildup and backscatter** effects
- Compare directly to **2D planar dose measurements** from QA devices

## Related Files

- **Configuration**: `tps_env/moqui_tps.in`
- **Implementation**: `moqui/base/environments/mqi_tps_env.hpp`
- **Grid geometry**: `moqui/base/mqi_grid3d.hpp`
- **Output handling**: `moqui/base/mqi_io.hpp`
