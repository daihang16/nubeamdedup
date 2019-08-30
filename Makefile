#***************************************************************************   
#       user configuration
# set mpi = yes for the parallel version
# set debug = yes to debug 
DEBUG ?= no  
# end user configuration;
#***************************************************************************    

#SHELL = /bin/sh

CC = g++

ifeq ($(strip $(DEBUG)), yes)
    CFLAGS += -ggdb -DDEBUG
else
    CFLAGS += -O2
endif

CFLAGS += -Wall -m64 -std=c++11 

STATICLIBS += /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a 

CPPFLAGS = CFLAGS

LIBS += -lm -lz #-lgsl -lgslcblas  

OBJS = model.o 

SRCS = model.cpp 


all:  model $(OBJS)

.cpp.o : ; 	$(CC) $(CFLAGS) -c $<  			

## static link to gsl libs.  		
## note it's crutial to use -static-libgcc instead of static; 
model: model.o $(OBJS); $(CC) -static-libgcc  $(CFLAGS) $(OBJS) $(STATICLIBS) $(LIBS) -o nubeamdedup	
#fp: fp.o $(OBJS); $(CC) -static-libgcc $(CFLAGS) $(OBJS) fp.o $(LIBS) -o bimbam	
#fp: fp.o $(OBJS); $(CC) $(CFLAGS) $(OBJS) fp.o $(STATICLIBS) $(LIBS) -o bimbam	

# the following dependence generated by make showdep	

showdep: 	
		@$(CC) -MM $(SRCS)
 
clean: 
		rm -f *.o
		
