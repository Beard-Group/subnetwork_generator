# Copyright 2009-2012 Karsten Ahnert
# Copyright 2009-2012 Mario Mulansky
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or
# copy at http://www.boost.org/LICENSE_1_0.txt)

#CUDA_ROOT = /usr/local/cuda-5.0
#ARCH = sm_35
HOME = /home/rthiagarajan
#BOOST_ROOT = /usr/include/boost
#THRUST_ROOT = /usr/local/cuda-5.0/include/thrust
#BOOST_ROOT = $(HOME)/include/boost/boost
SRC = .

#GCC = gcc-4.6.3
#CXX = gcc-4.6.3
GCC = gcc
CXX = gcc
#GCC = /usr/local/gcc-4.6.3/bin/g++-4.6.3

#INCLUDES += -I../../../../.. -I$(THRUST_ROOT) 

GCCFLAGS = -O3 -std=c++11 #$(INCLUDES)-Xcompiler -fopenmp -DTHRUST_HOST_BACKEND=THRUST_HOST_BACKEND_OMP
#NVCCFLAGS = -O3 $(INCLUDES) -arch $(ARCH) --compiler-bindir=/usr/local/gcc-4.6.3/bin/g++-4.6.3 -Xcompiler -DTHRUST_HOST_BACKEND=THRUST_HOST_BACKEND_CPP
#NVCCFLAGS = -O3 $(INCLUDES) --compiler-bindir=/home/raghu/gcc-4.6/bin/g++-4.6  -Xcompiler -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP

#LDLIBS = -lcudart -lgomp -lm -lstdc++
LDLIBS = -lm -lstdc++ 
#LDFLAGS = -L$(CUDA_ROOT)/lib64

%.co : %.cpp
	$(GCC) $(GCCFLAGS) -c $< -o $@

all : run_g645_250b

run_g645_250b:   $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/optimization/opt.co $(SRC)/subnetsearch/search.co main.co
#lsoda_run :   $(SRC)/io/io.co $(SRC)/lsoda/cuLsoda.co main.co
#lsoda_run :   $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co main.co
#run :   $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co $(SRC)/MCsim/mc.co $(SRC)/neldermead/nm.co main.co
	$(GCC) -o run_g645_250b $(LDLIBS) $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/optimization/opt.co $(SRC)/subnetsearch/search.co main.co
#	$(CC) -o lsoda_run $(LDFLAGS) $(LDLIBS) $(SRC)/lsoda/cuLsoda.co $(SRC)/io/io.co main.co
#	$(CC) -o run $(LDFLAGS) $(LDLIBS) $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co main.co
#	$(CC) -o run $(LDFLAGS) $(LDLIBS) $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co $(SRC)/MCsim/mc.co $(SRC)/neldermead/nm.co main.co
$(SRC)/io/io.co : $(SRC)/io/io.cpp
#$(SRC)/lsoda/cuLsoda.co : $(SRC)/lsoda/cuLsoda.cu
$(SRC)/spline/spline.co : $(SRC)/spline/spline.cpp
$(SRC)/optimization/opt.co : $(SRC)/optimization/opt.cpp
$(SRC)/subnetsearch/search.co : $(SRC)/subnetsearch/search.cpp
#$(SRC)/MCsim/mc.co : $(SRC)/MCsim/mc.cu
#$(SRC)/neldermead/nm.co : $(SRC)/neldermead/nm.cu
main.co : main.cpp

clean :
#	-rm *~ *.o *.co run
	rm -f  *.o  *.co  *~ $(SRC)/*/*.o $(SRC)/*/*.co  
