# Optimize
OPTIM = -O2

# The system defines used in executing make
#FF = f77 -C
FF = gfortran 

# Glueing all the things together
FFLAGS   = $(OPTIM) 
#-Wall
LFLAGS   = $(OPTIM) 
#-Wall


#
# Names
#
COMMONS        = configure.h
SRCMAKEDUST    = makedust.F 
OBJMAKEDUST    = makedust.o
SRCQMIE        = q_mie.F 
OBJQMIE        = q_mie.o
SRCQCDE        = q_cde.F 
OBJQCDE        = q_cde.o

OBJ            = $(OBJMAKEDUST) $(OBJQMIE) $(OBJQCDE)


#################################################
#                   RULES                       #
#################################################

all:	      makedust

makedust:     $(OBJ) Makefile
	      $(FF) $(LFLAGS) $(OBJ) $(LIBS) -o $@ 

makedust.o:   $(SRCMAKEDUST) Makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  makedust.F -o $@

q_mie.o:      $(SRCQMIE) Makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  q_mie.F -o $@

q_cde.o:      $(SRCQCDE) Makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  q_cde.F -o $@

clean:
	@rm -f	*.o *.mod
	@echo OBJECT and MODULE files removed.

cleanall:
	@rm -f	*.o *.mod *.dat *.out *.info *.used                 \
	   frequency.inp dustkappa_*.inp dust.inp                   \
	   makedust *~ *.lnk
	@echo Directory cleaned to basic
