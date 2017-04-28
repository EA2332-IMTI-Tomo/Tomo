## -*- makefile -*- 

# attention: -DNONE est utilisé par openCV et ca déconne

#GPU_SUPPORT = CUDA
GPU_SUPPORT = SOFT
SSE_SUPPORT = -DSSE4 -msse4  
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
OPTIMIZE =  -O3 #-pg
#-O3 $(MED_OPT) $(HARD_OPT) $(HYST_OPT) $(FINAL_OPT)


# si monocore
#FFTW = -lfftw3
# si multicore
FFTW_LIBS = -lfftw3_threads -lfftw3 -lfftw3f -lfftw3f_threads -lm -lpthread 
BOOST_LIBS = -lboost_system -lboost_chrono -lboost_thread -lboost_filesystem
OPENCV_LIBS = -lopencv_core -lopencv_highgui -lopencv_imgproc

EXTRA_LIBS = -lCamera -lPGM


COMMONFLAGS = -D$(GPU_SUPPORT) -D$(RELEASE_MODE) $(SSE_SUPPORT)  $(OPTIMIZE) $(INCL) $(DEBUG) 
#$(DEBUG)
SOFTFLAGS = -ansi -Wall -fopenmp #-fPIC
#-msse
# pour empêcher conversion silencieuse de double vers float

CFLAGS = $(COMMONFLAGS) $(SOFTFLAGS)



TOMO_LINK = $(FFTW_LIBS) $(BOOST_LIBS) $(OPENCV_LIBS) -lAIRVOL -lvectra  -lmsleep -fopenmp 
SHOW_LINK = $(TOMO_LINK) $(EXTRA_LIBS)

TOMO_OBJ_DIR = ../obj
TOMO_OBJS = $(TOMO_OBJ_DIR)/cpu_launcher.o $(TOMO_OBJ_DIR)/recon_parser.o $(TOMO_OBJ_DIR)/util.o $(TOMO_OBJ_DIR)/util_Image.o $(TOMO_OBJ_DIR)/Holo_Process.o $(TOMO_OBJ_DIR)/Compute.o $(TOMO_OBJ_DIR)/Focus_Corrector.o 

OWN_OWBS = $(TOMO_OBJS) 



OBJS = $(TOMO_OBJS) 

EXE = show_holo


#*******************************************************************************
# main section
#*******************************************************************************




all: $(OWN_OBJS) $(EXE)

$(EXE): $(OBJS) show_holo.o #*.h
	$(CC) $(OBJS) show_holo.o -o  $@ $(SHOW_LINK) 


#$^: toutes les dépendances
#$<: la première seulement


%.o: %.cc
	$(CC) -c $(CFLAGS)  -o $@ $< 


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
	rm -f *~  *.o $(EXE) Code.ps*

cleanobj:
	rm -f *~ $(OBJS) 

remake:
	make clean; make


tags:
	@etags  *.cc *.cu *.h	


.PHONY:all print display pdf clean remake tags
