CXX            = g++

CFLAGS        = -O2

FC            = gfortran

FFLAGS        = -O2

DEST	      = C:\Users\Tadaaki\Dropbox\Research\Libraries_64
#DEST         = C:\Users\Tadaaki\Dropbox\Research\SAT_2D\sat2d_ok\Oristano_example

INC_DIRS	  = C:\Users\Tadaaki\Dropbox\Research\Libraries_64\Headers

#HDRS	      = red_model.h addon.h

LDFLAGS	      = -L C:\Users\Tadaaki\Dropbox\Research\Libraries_64 

LIBS	      = -l blas -l lapack -l gfortran

LINKER	      = $(CXX)

MAKEFILE      = Makefile

#OBJS	      = axbsym.o pmat.o read.o red_model.o redmat.o svd.o
OBJS	      = matrix.o

PRINT	      = pr

CAT	          = cat

#PROGRAM	      = reduce_model.exe
PROGRAM	      = matrix

#SRCS	      = axbsym.cpp pmat.cpp read.cpp red_model.cpp redmat.cpp svd.cpp
#SRCS	      = build_red_model.cpp write.cpp read.cpp svd2.cpp

.SUFFIXES: .cpp .f

.cpp.o:
	$(CXX) $(CXXFLAGS) -I$(INC_DIRS) $(LDFLAGS) $(LIBS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

PRJS= $(SRCS:jo=.prj)

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS)
		ar rvs  $(DEST)\$(PROGRAM).lib $(OBJS)
#		$(CXX) $(OBJS)  -o $(PROGRAM).exe $(LDFLAGS) $(LIBS)

check: $(PRJS)
	ftnchek -noverbose -declare $(PRJS) -project -noextern -library > $(PROGRAM).ftn

profile:;       @make "CLAGS=-pg" "CFLAGS=-pg" "LDFLAGS=-pg" $(PROGRAM)

debug:;         @make "CFLAGS=-g" $(PROGRAM)

clean:;		rm *.o $(PROGRAM).lib

depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		$(DEST)\$(PROGRAM).lib 

print:;		$(PRINT) $(HDRS) $(SRCS)

cat:;		$(CAT) $(HDRS) $(SRCS) > $(PROGRAM).cat

program:        $(PROGRAM)

profile:        $(PROFILE)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

###
datin.o: SAT2D.H
ini2d.o: SAT2D.H
sat2d.o: SAT2D.H
satf2d.o: SAT2D.H
shl2d.o: SAT2D.H
