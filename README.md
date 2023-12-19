# Curvel (Curvelet) Construction in GPU
### Research @ LEMS, Brown University

This is a curvel/curvelet construction CPU and GPU implementation code. A "curvel" is a small segement of a curve constructed by a number of local third-order edges, capturing tangents, curvatures, etc. of the curve. The idea comes from the paper [Differential Geometry in Edge Detection](https://ieeexplore.ieee.org/abstract/document/8382271). For the third-order edge detection, refer to my [another GitHub](https://github.com/C-H-Chien/Third-Order-Edge-Detector) for more information. <br /> <br />

For CPU implementation, it was rewritten from the [official GitHub](https://github.com/yuliangguo/Differential_Geometry_in_Edge_Detection) of the paper, which basically replaced the data structure to organize the data and control the algorithm flow in order to copy necessary data from host (CPU) to device (GPU). The rewritten version has around 2 to 3 times faster than the original version, from 31,042 third-order edges in a 480x320 resolution of an image.

## Change Logs
June 15th, 2022: Complete rewriting the CPU version.
July 20th, 2022: First version of GPU implementation.

## Dependencies
The code is tested under a Linux-based system (specifically red-hat and Ubuntu). It might also work on other systems but I have not yet tested.
(1) cmake 3.15.4 or higher <br />
(2) gcc 10.2 or higher <br />
(3) cuda 11.1.1 or higher <br />

## How to use the code
(1) Compile the code simply by <br />
```bash
$ make curvelet
```
(2) Once the compilation is successful, you shall see an executable file which can be run by
```bash
$ ./CURVELET
```

