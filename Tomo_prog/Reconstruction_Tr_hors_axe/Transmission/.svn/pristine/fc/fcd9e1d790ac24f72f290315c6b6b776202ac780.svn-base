## -*- makefile -*- 

# attention: -DNONE est utilisé par openCV et ca déconne

#GPU_SUPPORT = CUDA
GPU_SUPPORT = SOFT
RELEASE_MODE = DEBUG_MODE
#RELEASE_MODE = FINAL

#NVCC = nvcc
CC = g++



# good optimization flags
MED_OPT = -finline-functions -fcaller-saves  -fthread-jumps -fstrength-reduce
# expensive optimizations
HARD_OPT = -funroll-loops -fexpensive-optimizations -frerun-cse-after-loop
# hysterical /#fastmath: loss of IEEE complying #-felide-constructors: can alter side effects of constructors
HYST_OPT = -ffast-math -felide-constructors
# omitfp: debugging impossible
FINAL_OPT = -fomit-frame-pointer 


DEBUG = -g
#OPTIMIZE = -g 
OPTIMIZE = -O3 -pg 
#-O3 $(MED_OPT) $(HARD_OPT) $(HYST_OPT) $(FINAL_OPT)


# si monocore
#FFTW = -lfftw3
# si multicore
FFTW_LIBS = -lfftw3_threads -lfftw3 -lfftw3f -lfftw3f_threads -lm -lpthread 
BOOST_LIBS = -lboost_system -lboost_chrono
OPENCV_LIBS = -lopencv_core -lopencv_highgui -lopencv_imgproc



COMMONFLAGS = -D$(GPU_SUPPORT) -D$(RELEASE_MODE) $(OPTIMIZE)  $(INCL)
#$(DEBUG)
SOFTFLAGS = -fPIC -ansi -Wall -fopenmp 
#-msse
# pour empêcher conversion silencieuse de double vers float
CUDAFLAGS = #--gpu-name sm_13
# --compiler-options -fpermissive // un peu brutos encore
#-use_fast_math -Xcompiler -D__builtin_stdarg_start=__builtin_va_start
NVCCFLAGS = $(COMMONFLAGS) $(CUDAFLAGS)
CFLAGS = $(COMMONFLAGS) $(SOFTFLAGS)


CUDA_LINK = -lcufft -lcuda 

TOMO_LINK = $(FFTW_LIBS) $(BOOST_LIBS) $(OPENCV_LIBS) -lAIRVOL -lvectra -fopenmp



SOURCES = fonctions.cpp main.cpp

OWN_OBJS = fonctions.o main.o

# $AIRVOL_OBJ   



OBJS = $(OWN_OBJS) 

EXE = transmission_soft
#EXE = tomoTrans_CUDA


#*******************************************************************************
# main section
#*******************************************************************************




all: $(OWN_OBJS) $(EXE)

$(EXE): $(OBJS) *.h
	$(CC) $(OBJS) -o  $@ $(TOMO_LINK) -pg


#$^: toutes les dépendances
#$<: la première seulement


%.o: %.cc
	$(CC) -c $(CFLAGS)  -o $@ $< 


#cuFourier.o: cuFourier.cu *.h
#	$(CC) $(CFLAGS) -c cuFourier.cu 




#*******************************************************************************
# phony utils
#*******************************************************************************


# Code printout/display
print:	
	a2ps -2  \
		main.h main.cc \
		-o Code.ps
	gzip Code.ps


display: Code.ps.gz
	gv Code.ps.gz &	


pdf: 	print
	gunzip Code.ps.gz
	ps2pdf Code.ps Code.pdf
# -dCompatibilityLevel=1.3 -dSubsetFonts=true -dEmbedAllFonts=true -sPAPERSIZE=letter
	gzip Code.ps


clean:
	rm -f *~ $(OBJS) $(EXE) Code.ps*


remake:
	make clean; make


tags:
	@etags  *.cc *.cu *.h	


.PHONY:all print display pdf clean remake tags
