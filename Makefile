CXX = g++
BLAS_INC = /home/mathlib/OpenBLAS/build/include
BLAS_LIB = /home/mathlib/OpenBLAS/build/lib/ -lopenblas
LIBS = -L $(BLAS_LIB)
LIBS += -ldl -lstdc++
CXXFLAGS = -I$(BLAS_INC) $(LIBS)

cp2k-tester : cp2k-tester.cpp
	$(CXX) -o cp2k-tester cp2k-tester.cpp $(CXXFLAGS)

clean : 
	-rm cp2k-tester
