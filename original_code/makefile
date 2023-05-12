#CUDADIR ?= /gpfs/runtime/opt/cuda/11.1.1/cuda
#INC=-I${CUDADIR}/include
#LIBDIR=-L${CUDADIR}/lib64
#LIB=-lcudart

#OPENMP = -fopenmp

#NVCCFLAGS  = -O3  -Xcompiler "-fPIC -Wall -Wno-unused-function -Wno-strict-aliasing" -std=c++11
#NVCCFLAGS += -gencode arch=compute_70,code=sm_70
#NVCCFLAGS += -gencode arch=compute_80,code=sm_80

curvelet: main.o CC_curve_model_3d.o curvelet.o curveletmap.o form_curvelet_process.o  
	g++ -o CURVELET main.o CC_curve_model_3d.o curvelet.o curveletmap.o form_curvelet_process.o  

main.o: main.cpp
	g++ -c main.cpp -O3 -o main.o

CC_curve_model_3d.o: CC_curve_model_3d.cpp
	g++ ${OPENMP} -c CC_curve_model_3d.cpp -O3 -o CC_curve_model_3d.o

curveletmap.o: curveletmap.cpp
	g++ ${OPENMP} -c curveletmap.cpp -O3 -o curveletmap.o

curvelet.o: curvelet.cpp
	g++ ${OPENMP} -c curvelet.cpp -O3 -o curvelet.o

form_curvelet_process.o: form_curvelet_process.cpp
	g++ ${OPENMP} -c form_curvelet_process.cpp -O3 -o form_curvelet_process.o

clean:
	rm -f curvelet *.o


