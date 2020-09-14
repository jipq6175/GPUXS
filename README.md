# GPUSWAXS
Calculating SWAXS profiles from large number of atoms using GPU

(Last Updated: 09/02/2019, Merged to SWAXS.jl, Archived 09/14/2020)


## Prerequisite
1. ArrayFire
2. Microsoft Visual Studio Community

## Usage

Run `gpuswaxs --help`

```

********************************************************************************
**********     gpuswaxs.exe: version 1.0.0                            **********
**********     A GPU-accelerated computation of the small-wide-angle  **********
**********     X-ray scattering (SWAXS) profile (.dat) from .binvox   **********
**********     or .pdb file(s). Require: ArrayFire Library            **********
**********                                                            **********
**********     AUTHOR: Yen-Lin Chen @ Cornell                         **********
**********     DATE:   Sept. 2 2019                                   **********
********************************************************************************

GPU implementation of SWAXS profile calculation.

Usage:
  GPUSWAXS [OPTION...] [ARGS...]

  -B, --binvox arg    The shape file of dummy voxels: .binvox (default: )
  -P, --pdb arg       The single PDB mode for atomic coordinates: .pdb (default: )
  -T, --solute arg    The solute.pdb file (default: )
  -V, --solvent arg   The solvent.pdb file (default: )
  -d, --density arg   The electron density averaged on a voxel, in e/A^3 (default: 0.50)
  -v, --voxdim arg    The voxel size in A (default: 2.0)
  -J, arg             The number of orientations (default: 1800)
      --qmin arg      The minimum q in 1/A
      --qspacing arg  The spacing in q axis
      --qmax arg      The maximum q reached
  -h, --help          Print this help message and exit

```


## Example

### 1. Calculating SWAXS profiles from a shape `.binvox`

   Run the following

```
gpuswaxs --binvox por.binvox 0.0 0.01 1.0
```

   The command line output should look like this:

```

********************************************************************************
**********     gpuswaxs.exe: version 1.0.0                            **********
**********     A GPU-accelerated computation of the small-wide-angle  **********
**********     X-ray scattering (SWAXS) profile (.dat) from .binvox   **********
**********     or .pdb file(s). Require: ArrayFire Library            **********
**********                                                            **********
**********     AUTHOR: Yen-Lin Chen @ Cornell                         **********
**********     DATE:   Sept. 2 2019                                   **********
********************************************************************************

GPUSWAXS: Operating in binvox shape mode ...
GPUSWAXS: binvox file = por.binvox with J = 1800 and q = 0:0.01:1 density = 0.5 voxdim = 2
Voxel INFO: [#binvox]
Voxel INFO: Reading binvox version: 1
Voxel INFO: Keyword [translate], ignored.
Voxel INFO: Keyword [scale], ignored.
Voxel INFO: Success, read 4083 voxels
Voxel INFO: por.binvox read-in successfully.
Voxel INFO: Start to calculate swaxs curves with J = 1800 and q = 0:0.01:1 density = 4e/A^3
100% completed: |||||||||||||||||||||||||||||||||||||||||||||||||
Voxel INFO: Elapsed time: 3.95272 seconds
Voxel INFO: Writing q, and intensity to the file: por.dat .
Voxel INFO: Writing completed successfully.

```

   The SWAXS profile is then saved as `por.dat` file.



### 2. Calculating SWAXS profile from single `.pdb` file of atomic coordinates

   Run the following

```
gpuswaxs --pdb rna.pdb 0.0 0.01 1.0
```

   The output should look like

```

********************************************************************************
**********     gpuswaxs.exe: version 1.0.0                            **********
**********     A GPU-accelerated computation of the small-wide-angle  **********
**********     X-ray scattering (SWAXS) profile (.dat) from .binvox   **********
**********     or .pdb file(s). Require: ArrayFire Library            **********
**********                                                            **********
**********     AUTHOR: Yen-Lin Chen @ Cornell                         **********
**********     DATE:   Sept. 2 2019                                   **********
********************************************************************************

GPUSWAXS: Operating in single pdb mode ...
GPUSWAXS: pdb file = rna.pdb with J = 1800 and q = 0:0.01:1
PDB INFO: Detected 8998 atoms in the file: rna.pdb.
PDB INFO: Processing all atomic coordinates ...
PDB INFO: Finished processing all atomic coordinates.
PDB INFO: Start to calculate swaxs curves...
100% completed: |||||||||||||||||||||||||||||||||||||||||||||||||
PDB INFO: Elapsed time: 8.91705 seconds
PDB INFO: Writing q, and intensity profile to the file: rna.dat .
PDB INFO: Writing completed successfully.

```

   The SWAXS profile is saved as `rna.dat` file.



### 3. Calculating SWAXS profile from solute and solvent pair

   Run

```
gpuswaxs --solute solute.pdb --solvent solvent.pdb 0.0 0.01 1.0
```

   The output is then

```

********************************************************************************
**********     gpuswaxs.exe: version 1.0.0                            **********
**********     A GPU-accelerated computation of the small-wide-angle  **********
**********     X-ray scattering (SWAXS) profile (.dat) from .binvox   **********
**********     or .pdb file(s). Require: ArrayFire Library            **********
**********                                                            **********
**********     AUTHOR: Yen-Lin Chen @ Cornell                         **********
**********     DATE:   Sept. 2 2019                                   **********
********************************************************************************

GPUSWAXS: Operating in pair pdb mode ...
GPUSWAXS: solute file = solute.pdb solvent file = solvent.pdb with J = 1800 and q = 0:0.01:1
PDB INFO: Detected 8998 atoms in the file: solute.pdb.
PDB INFO: Processing all atomic coordinates ...
PDB INFO: Finished processing all atomic coordinates.
PDB INFO: Detected 8230 atoms in the file: solvent.pdb.
PDB INFO: Processing all atomic coordinates ...
PDB INFO: Finished processing all atomic coordinates.
PDB INFO: Start to calculate swaxs curves...
100% completed: |||||||||||||||||||||||||||||||||||||||||||||||||
PDB INFO: Elapsed time: 13.8732 seconds
PDB INFO: Writing q, and intensity to the file: solute_buffer_subtracted.dat .
PDB INFO: Writing completed successfully.

```

   The result will have a postfix of `_buffer_subtracted`. The SWAXS file is saved as `solute_buffer_subtracted.dat`.



## Notes
1. Regarding PDB errors, please format your `.pdb` file(s).
2. Not sure if I will maintain this in the future ...
