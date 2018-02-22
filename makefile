# for compilation on Linux and Mac desktops:
# 1. use mpic++ compiler, and -fopenmp directive
CXX = mpic++
OMP = -fopenmp

# for compilation on CSU ISTEC Cray 
# 1. make clean
# 2. use CC compiler, and -fopenmp directive
#CXX = CC
#OMP = -fopenmp

# for compilation on TACC Lonestar
# 1. make clean
# 2. use mpicxx compiler, and -openmp directive
#CXX = mpicxx
#OMP = -openmp

FLAGS = -O3 -Wall

all: m+.o mp.o aStar.o 
	${CXX} ${FLAGS} -o m+ m+.o mp.o aStar.o ${OMP}

m+.o: m+.cpp m+.hpp
	${CXX} ${FLAGS} -c m+.cpp

mp.o: mp.cpp m+.hpp
	${CXX} ${FLAGS} -c mp.cpp

aStar.o: aStar.cpp m+.hpp
	${CXX} ${FLAGS} -c aStar.cpp ${OMP}

clean:
	rm -rf *.o
