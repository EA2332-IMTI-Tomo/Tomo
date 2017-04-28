_DEBUG = true



#c'est pas nécessaire car l'include joint le source
LIBVECTRA = lib/libvectra.so


EXECS     =  $(LIBVECTRA)


#ifeq ($(SRC_CS) $(SRC_CPPS),)
#  $(error No source files specified)
#endif

#ifeq ($(EXEC),)
#  $(error No executable file specified)
#endif

CC=gcc
CPP=g++
LD=$(CPP)



INC_DIR         = include
LIB_DIR         = lib




ifdef _DEBUG
CFLAGS    += -g -D_DEBUG
CPPFLAGS  += -g -D_DEBUG
else
CFLAGS    += -O3
CPPFLAGS  += -O3
endif


BOOST_LIBS = -lboost_system -lboost_chrono

CFLAGS    += -I$(INC_DIR) -D_UNIX_ -D_LINUX_
CPPFLAGS  += -I$(INC_DIR) -D_UNIX_ -D_LINUX_
LDFLAGS   += -L$(LIB_DIR)         \
	     $(BOOST_LIBS) 	  

OBJS      += $(SRC_CPPS:%.cc=%.o)
OBJS      += $(SRC_CS:%.c=%.o)

all: $(EXECS) 


$(LIBVECTRA): vectra.o
	$(CPP) $< -shared -o $@ $(LDFLAGS)

vectra.o: vectra.cc
	$(CPP) -c $< $(CPPFLAGS) -fPIC -o $@ 

#bin/manipTomoHDR: manipTomoHDR.o
#	$(CPP) $<  -o $@  $(LDFLAGS) $(OPENCV_LINK)



%.o: %.cc
	$(CPP) -c $(CPPFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

clean:
	rm -rf *.o $(EXECS)

tags: 
	@etags *.cc *.c *.h --lang=c++


.PHONY: all clean tags
