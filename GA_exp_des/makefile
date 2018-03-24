# -*- Mode: makefile -*-
# Copyright (c) 1999 Matthew Wall, all rights reserved
# -----------------------------------------------------------------------------
#   To make all of the examples, do 'make'.  You can compile any one of
# the examples by typing 'make exN' where N is the number of the example you
# want to compile.  See the README for a description of what each example does.
# -----------------------------------------------------------------------------

include C:\Users\Tadaaki\Dropbox\Research\galib247\makevars
#include C:\Users\Tadaaki\Dropbox\Research\galib247\makevars
FC				= gfortran
FFLAGS        	= -O2
# Set these paths to the location of the GA library and headers.
#GA_INC_DIR= /usr/local/include
#GA_LIB_DIR= /usr/local/lib

###32 bit libraries
#GA_INC_DIR		= C:\Users\Tadaaki\Dropbox\Research\galib247
#GA_LIB_DIR		= C:\Users\Tadaaki\Dropbox\Research\galib247\ga
#LAPACK_LIB_DIR	= C:\Users\Tadaaki\Dropbox\Research\Libraries
#HDR_DIRS		= C:\Users\Tadaaki\Dropbox\Research\Libraries\Headers

###64 bit libraries
GA_INC_DIR		= C:\Users\Tadaaki\Dropbox\Research\Libraries_64\galib247

GA_LIB_DIR		= C:\Users\Tadaaki\Dropbox\Research\Libraries_64\galib247\ga

LAPACK_LIB_DIR	= C:\Users\Tadaaki\Dropbox\Research\Libraries_64\

HDR_DIRS		= C:\Users\Tadaaki\Dropbox\Research\Libraries_64\Headers



GA_libs1 		= ga
GA_libs2 		= m

INC_DIRS		= -I$(HDR_DIRS) -I$(GA_INC_DIR)
LIB_DIRS		= -L$(GA_LIB_DIR) -L$(LAPACK_LIB_DIR)
libs 			= -l$(GA_libs1) -l$(GA_libs2) -l matrix -l blas -l lapack   -l gfortran

EXS 			= Exp_Des_GA_param
OBJS	      	= GAsearch.o read.o  axbsym.o write.o\
				  red_solver2.o red_mat.o obj.o dist2.o\
				  fortsys.o\

#OBJS			= test.o\
#				  read.o  axbsym.o write.o\
#				  red_solver2.o red_mat.o obj.o dist2.o fortsys.o\



#DEST			= C:\Users\Tadaaki\Dropbox\Research\SAT_2D\sat2d_ok\Oristano_example
#DEST	    	= C:\Users\Tadaaki\Dropbox\Research\SAT_2D\rectangular_example
DEST			= C:\Users\Tadaaki\Dropbox\Research\SAT_2D\Medium_model\med_rect

.SUFFIXES: .cpp .f
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC_DIRS) -std=c++11  -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(EXS)
# Use this for gnu make
#$(EXS): %: %.o
#	$(CXX) $@.o -o $@ $(LIB_DIRS) $(libs) $(CXX_LIBS)

$(EXS):  $(OBJS)
	$(CXX)  $(OBJS) -o $(EXS).exe $(LIB_DIRS) $(libs) $(CXX_LIBS)

clean:
	$(RM) $(EXS)
	$(RM) *.o *~ *.bak *.pixie core
	$(RM) test_results.txt test_stats.txt
	$(RM) *.dat
	$(RM) *.out *.exe vc* *.pdb

test: $(EXS)
	$(RM) test_results.txt test_stats.txt
	@echo "running tests.  this could take up to 1/2 hour, depending on"
	@echo "the speed of your computer.  monitor test_results.txt and"
	@echo "test_stats.txt to see what is happening."
	@echo ""
	@rm -f test_results.txt
	@echo `uname -a` > test_stats.txt
	@echo "" >> test_stats.txt
	for x in $(EXS); do \
          echo "$$x... "; \
          echo "$$x" >> test_stats.txt; \
          echo "start:  " `date` >> test_stats.txt; \
          ./$$x seed 555 >> test_results.txt; \
          echo "finish: " `date` >> test_stats.txt; \
          echo "" >> test_stats.txt; \
        done

install: $(EXS)
		install -s $(EXS) $(DEST)