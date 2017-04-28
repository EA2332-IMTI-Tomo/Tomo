CC = g++

BASE_DIR = /users/these/bailleul/These

AIR_VOL_DIR = $(BASE_DIR)/AIR_VOL
AIR_VOL_INCL = -I$(AIR_VOL_DIR)
AIR_VOL_OBJS = $(AIR_VOL_DIR)/Filename.o

OPENSOURCE_DIR = $(BASE_DIR)/OpenSource
OPENSOURCE_INCL = -I$(OPENSOURCE_DIR)
OPENSOURCE_OBJS = $(OPENSOURCE_DIR)/Tools.o

INCL = $(OPENSOURCE_INCL) $(AIR_VOL_INCL)
CFLAGS = -g -ansi -Wall $(INCL)

OWN_OBJS = 
OBJS = $(OWN_OBJS) $(AIR_VOL_OBJS) $(OPENSOURCE_OBJS)

EXE = opensource


$(EXE): $(OBJS) 

Tools.o: Tools.c
	$(CC) $(CFLAGS) -c Tools.c