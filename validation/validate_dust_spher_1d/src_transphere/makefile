# Optimize
#OPTIM = -O2

# The system defines used in executing make
#FF = f77 
FF = g77 

# Glueing all the things together
FFLAGS   = $(OPTIM) 
#-Wall
LFLAGS   = $(OPTIM) 
#-Wall


#
# Names
#
COMMONS        = common_transfer.h common_radtrans.h common_dust.h configure.h
SRCVETRADTRANS = vetradtrans.F 
OBJVETRADTRANS = vetradtrans.o
SRCRADMOMEQ    = radmomeq.F 
OBJRADMOMEQ    = radmomeq.o
SRCTRANSFER    = transfer.F 
OBJTRANSFER    = transfer.o
SRCDUST        = dust.F 
OBJDUST        = dust.o
SRCNRECIP      = nrecip.F 
OBJNRECIP      = nrecip.o
SRCMAIN        = main.F 
OBJMAIN        = main.o

OBJ            = $(OBJVETRADTRANS) $(OBJRADMOMEQ)  $(OBJTRANSFER)  $(OBJMAIN)  $(OBJDUST)  $(OBJNRECIP) 


#################################################
#                   RULES                       #
#################################################

all:	      transphere

transphere:   $(OBJ) makefile
	      $(FF) $(LFLAGS) $(OBJ) $(LIBS) -o $@ 

vetradtrans.o:$(SRCVETRADTRANS) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  vetradtrans.F -o $@

radmomeq.o:   $(SRCRADMOMEQ) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  radmomeq.F -o $@

transfer.o:   $(SRCTRANSFER) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  transfer.F -o $@

dust.o:       $(SRCDUST) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  dust.F -o $@

nrecip.o:     $(SRCNRECIP) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  nrecip.F -o $@

main.o:       $(SRCMAIN) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  main.F -o $@
