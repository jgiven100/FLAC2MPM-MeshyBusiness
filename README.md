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
Example input file:
```
{
    "xmin": -100,
    "xmax": 100,
    "ymin": 0,
    "ymax": 50,
    "he": 1.0,
    "ppc": 2,
    "file": "FLAC-output.f2grid"
}
```

Definitions:
* `"xmin"` : Minimum x-value for mesh
* `"xmax"` : Maximum x-value for mesh
* `"ymin"` : Minimum y-value for mesh
* `"ymax"` : Maximum y-value for mesh
* `"he"` : Characteristic length (height) of element
* `"ppc"` : Particles per cell where `PPC=N` means `N x N` distribution (equally spaced)
* `"file"` : Name of FLAC output file

## Outputs

* Mesh file
* Entity sets json
    * `nset_id=0` := left boundary
    * `nset_id=1` := right boundary
    * `nset_id=2` := bottom boundary
* Particle coordinates
* Particle cells
* Particle volumes

*TODO*
* Particle stresses
* Particle stresses (beginning)
