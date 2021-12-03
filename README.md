# Framework for surface mesh processing

The program is used to load meshes and render using OpenGL.

## External Libraries

* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/), Recommended version: the latest 8.1(at Oct. 2020)
* [Qt](https://www.qt.io/), Recommended version: 5.13.0

## Usage

```
git clone https://github.com/USTC-GCL-F/Surface-Mesh-Framework
cd SurfaceFramework
```

Edit lines 7-9 of CmakeLists.txt to set the values of **EIGEN_PATH**,**OPENMESH_PATH** and **OPENMESH_LIB_PATH**
```
mkdir build && cd build
cmake -A x64 ..
```

Open **SurfaceFramework.sln**, select **SrufaceFramework** as launch project, and run.


## Supported File Formats

.obj .off .ply .stl

## QWQ

这里存放了DGP、FEM和CAGD部分作业，还有一些DGP的作业下次一定挪过来qwq

目前包含

ARAP的曲面形变->APAPSurfaceModeling.h

MVC参数化->MVAParameterization.h

2d有限元求解pde->FEM2D.h

Loop细分->LoopSubdivision.h