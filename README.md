# Curvel (Curvelet) Construction in GPU
### Research @ LEMS, Brown University

This is a curvel/curvelet construction CPU and GPU implementation code. A "curvel" is a small segement of a curve constructed by a number of local third-order edges, capturing tangents, curvatures, etc. of the curve. The idea comes from the paper [Differential Geometry in Edge Detection](https://ieeexplore.ieee.org/abstract/document/8382271). For the third-order edge detection, refer to my [another GitHub](https://github.com/C-H-Chien/Third-Order-Edge-Detector) for more information. <br />

For CPU implementation (``cpu_curvelet`` folder, with optional OpenMP parallelism via ``--nthreads N``), it was rewritten from the code (see ``original_code`` directory) released in the [official GitHub](https://github.com/yuliangguo/Differential_Geometry_in_Edge_Detection) of the paper. The rewritten version was basically replacing the data structure to organize the data and control the algorithm flow. 

## Repository layout

| Folder | Description |
|--------|------|
| `original_code/` | Reference C++ from the official released code of the paper |
| `cpu_curvelet/` | Rewritten CPU version (flat arrays for edges and curve bundles) with OpenMP support |
| `gpu_curvelet/` | (TBD) GPU port of the rewrite |
| `test_files/` | Sample third-order edge files |

The rewrite (`cpu_curvelet`) keeps the same core math and replaces object-heavy curve-bundle and curvelet storage (`edgemap`, `curveletmap`, `CC_curve_model_3d`) with flat arrays in the construction kernel. This was intended for the friendliness of GPU implementation.

## Dependencies
The code is tested under a Linux-based system (specifically red-hat and Ubuntu). It might also work on other systems but I have not yet tested. <br />
(1) cmake 3.15.4 or higher <br />
(2) gcc 10.2 or higher <br />
(3) cuda 11.1.1 or higher for GPU impelemtation <br />

## How to use the code
(1) Navigate to either ``cpu_curvelet`` or `original_code`. <br />
(2) Compile the code simply by <br />
```bash
$ make curvelet
```
(2) Once the compilation is successful, you shall see an executable file which can be run by
```bash
$ ./CURVELET
```
### Inputs
There are multiple input arguments allowed to be specified at runtime. Their default values are defined in `cpu_curvelet/param_settings.hpp` or `gpu_curvelet/param_settings.hpp`. The mostly used input argument is the third-order edges file:
```bash
$ ./CURVELET --edge-file your-edge-file.txt
```
Example edge files live under `test_files/`. For CPU implementation, OpenMP parallelism is supported. The number of cores is 1 by default, which can be changed by setting `nthreads`, _i.e._,
```bash
$ ./CURVELET --nthreads 8
```

### Outputs
Output chains and curvelet information are written under `outputs/`.

## Workflow of Curvelet Construction
See the documentation of the high-level overview of the curvelet construction [here](https://github.com/C-H-Chien/Curvelet_Construction_on_GPU/blob/init_dev/docs/curvelet_construction.md).

## TODO
More updates and documentations are coming soon.
