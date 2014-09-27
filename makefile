PROG		= main.exe
#OBJCTS		= globvar.o func.o sub.o solver.o main.o
OBJCTS		= func.o sub.o solver.o main.o

FC		= gcc

## normal ver.
FFLAGS		= -O4 #-ffloat-store

## debug ver.
#FFLAGS		= -CB -traceback -g -heap-arrays #-check uninit -warn all -check all -std
#FFLAGS		= -O0 -g -Wall # -fbounds-check -fbacktrace -O -Wuninitialized

#COMMON_MOD 	= globvar.f90 func.f90 solver.f90
COMMON_MOD 	= func.c sub.c solver.c

.SUFFIXES	:
.SUFFIXES	: .c .o

.c.o:
	${FC} -c $< ${FFLAGS} 

${PROG}	:	${OBJCTS}
	${FC} -o $@ ${OBJCTS} ${FFLAGS} >& err.d

${OBJCTS}: ${COMMON_MOD}

clean:
	rm -f ${PROG} *.o err.d info.dat reshis.dat solution.dat
