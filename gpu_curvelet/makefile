CUDADIR ?= $(if $(CUDA_HOME),$(CUDA_HOME),/gpfs/runtime/opt/cuda/11.1.1/cuda)
INC=-I${CUDADIR}/include
LIBDIR=-L${CUDADIR}/lib64
LIB=-lcudart

#OPENMP = -fopenmp

NVCCFLAGS  = -O3  -Xcompiler "-fPIC -Wall -Wno-unused-function -Wno-strict-aliasing" -std=c++11
NVCCFLAGS += -gencode arch=compute_70,code=sm_70
NVCCFLAGS += -gencode arch=compute_80,code=sm_80

CXXFLAGS = -O3 -std=c++17

curvelet: main.o gpu_curvelet.o gpu_preprocess.o
	g++ ${INC} ${CXXFLAGS} -o CURVELET main.o gpu_curvelet.o gpu_preprocess.o ${LIBDIR} ${LIB}

main.o: main.cpp data_io.hpp param_settings.hpp
	g++ ${INC} ${CXXFLAGS} -c main.cpp -o main.o

gpu_curvelet.o: gpu_build_curvelet.cu
	nvcc ${INC} ${NVCCFLAGS} -c gpu_build_curvelet.cu -o gpu_curvelet.o

gpu_preprocess.o: gpu_preprocess.cu gpu_preprocess.hpp gpu_dev_preprocess.cuh gpu_kernels.hpp
	nvcc ${INC} ${NVCCFLAGS} -c gpu_preprocess.cu -o gpu_preprocess.o

clean:
	rm -f curvelet *.o


