# FLAC2MPM-MeshyBusiness
Generate mesh and particle input files for MPM runout analysis based on FLAC output.

## Compile
Start from `FLAC2MPM-MeshyBusiness` directory:
```
mkdir build/
g++ main.cpp -std=c++17 -Wall -o build/main
```


## Run
Copy `input.json` and `FLAC-output.f2grid` files in the build directory:
```
cd build/
./main -i input.json
```

## Input file
**Example input file:**
```
{
    "xmin": -100,
    "xmax": 100,
    "ymin": 0,
    "ymax": 50,
    "he": 1.0,
    "ppc": 2,
    "file": "FLAC-output.f2grid"
    "stress_file": "FLAC-output-stress.txt"
    "material_file": "FLAC-output-material.txt"
}
```

**Definitions:**
* `"xmin"` : Minimum x-value for mesh
* `"xmax"` : Maximum x-value for mesh
* `"ymin"` : Minimum y-value for mesh
* `"ymax"` : Maximum y-value for mesh
* `"he"` : Characteristic length (height) of element
* `"ppc"` : Particles per cell where `PPC=N` means `N x N` distribution (equally spaced)
* `"file"` : Name of FLAC output file
* `"stress_file"` : Name of FLAC output stress file

**Stress file:**

Stress file should be a text file comprised of rows with zone id, x-dir normal stress (effective), y-dir normal stress (effective), xy-dir shear stress, and pore pressure. Do not include a header. Use `space` as a delimiter. Stresses and pressure should use *tension positive*.
```
zid  sigma'_xx  sigma'_yy  sigma_xy  u
```

**Material file:**

Material file should be a text file comprised of rows with zone id, bulk modulus, shear modulus, density, liquefied strength, and drained strength. Do not include a header. Use `space` as a delimiter. Drained strength is adopted unless liquefied strength `!=0`.
```
zid  shear_modulus  bulk_modulus  density  c_liq  c_dry
```

## Outputs

* Mesh file
* Entity sets json
    * `nset_id=0` := left boundary
    * `nset_id=1` := right boundary
    * `nset_id=2` := bottom boundary
* Particle coordinates
* Particle cells
* Particle volumes
* Particle stresses
* Particle stresses (effective)
* Particle materials

## Python Preproc

The "preproc.py" code converts the FLAC output .csv files into the stress.txt and material.txt files which are inputs to the cpp code described above.
