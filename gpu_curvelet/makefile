CUDADIR ?= $(if $(CUDA_HOME),$(CUDA_HOME),/gpfs/runtime/opt/cuda/11.1.1/cuda)
INC=-I${CUDADIR}/include -I.
LIBDIR=-L${CUDADIR}/lib64
LIB=-lcudart

NVCCFLAGS  = -O3 -Xcompiler "-fPIC -Wall -Wno-unused-function -Wno-strict-aliasing" -std=c++17
NVCCFLAGS += -gencode arch=compute_70,code=sm_70
NVCCFLAGS += -gencode arch=compute_80,code=sm_80

CXXFLAGS = -O3 -std=c++17 -I.

curvelet: main.o gpu_preprocess.o gpu_curve_bundle_formation.o
	g++ ${INC} ${CXXFLAGS} -o CURVELET main.o gpu_preprocess.o gpu_curve_bundle_formation.o ${LIBDIR} ${LIB}

main.o: main.cpp data_io.hpp param_settings.hpp gpu_common.hpp gpu_curve_bundle_formation.hpp
	g++ ${INC} ${CXXFLAGS} -c main.cpp -o main.o

gpu_curve_bundle_formation.o: gpu_curve_bundle_formation.cu gpu_curve_bundle_formation.hpp gpu_common.hpp timer.hpp \
		device/gpu_curve_bundle_formation_kernels.cuh gpu_dev_curve_bundle.cuh indices.hpp
	nvcc ${INC} ${NVCCFLAGS} -c gpu_curve_bundle_formation.cu -o gpu_curve_bundle_formation.o

gpu_preprocess.o: gpu_preprocess.cu gpu_preprocess.hpp gpu_common.hpp \
		device/gpu_dev_common.cuh device/gpu_dev_neighbors.cuh device/gpu_preprocess_kernels.cuh
	nvcc ${INC} ${NVCCFLAGS} -c gpu_preprocess.cu -o gpu_preprocess.o

clean:
	rm -f CURVELET *.o
