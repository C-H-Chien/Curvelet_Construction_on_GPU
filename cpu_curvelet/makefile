#CUDADIR ?= /gpfs/runtime/opt/cuda/11.1.1/cuda
#INC=-I${CUDADIR}/include
#LIBDIR=-L${CUDADIR}/lib64
#LIB=-lcudart

#OPENMP = -fopenmp

#NVCCFLAGS  = -O3  -Xcompiler "-fPIC -Wall -Wno-unused-function -Wno-strict-aliasing" -std=c++11
#NVCCFLAGS += -gencode arch=compute_70,code=sm_70
#NVCCFLAGS += -gencode arch=compute_80,code=sm_80

curvelet: main.o   
	g++ -o CURVELET main.o 

main.o: main.cpp
	g++ -c main.cpp -O3 -o main.o

clean:
	rm -f curvelet *.o


