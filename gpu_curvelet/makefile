CUDADIR ?= /gpfs/runtime/opt/cuda/11.1.1/cuda
INC=-I${CUDADIR}/include
LIBDIR=-L${CUDADIR}/lib64
LIB=-lcudart

#OPENMP = -fopenmp

NVCCFLAGS  = -O3  -Xcompiler "-fPIC -Wall -Wno-unused-function -Wno-strict-aliasing" -std=c++11
NVCCFLAGS += -gencode arch=compute_70,code=sm_70
NVCCFLAGS += -gencode arch=compute_80,code=sm_80

curvelet: main.o gpu_curvelet.o   
	g++ ${INC} -o CURVELET main.o gpu_curvelet.o ${LIBDIR} ${LIB}

main.o: main.cpp
	g++ -c main.cpp -O3 -o main.o

gpu_curvelet.o: gpu_build_curvelet.cu
	nvcc ${INC} ${NVCCFLAGS} -c gpu_build_curvelet.cu -o gpu_curvelet.o

clean:
	rm -f curvelet *.o


