BestFitGeoplanes.out: InvertRoutine.o InvertGeoplanes.o
	gfortran InvertRoutine.o InvertGeoplanes.o -lstdc++ -o BestFitGeoplanes.out -L/usr/lib/lapack -llapack  -L/usr/lib/libblas -lblas 

InvertRoutine.o bestfitplane.mod geologic_processing.mod math_constants.mod var_types.mod vector_processing.mod: InvertRoutine.f03
	gfortran -c InvertRoutine.f03 

InvertGeoplanes.o: InvertGeoplanes.cpp
	g++ -c -std=c++11 InvertGeoplanes.cpp

.PHONY: clean
clean:
	rm *.mod *.o


