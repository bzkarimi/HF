FC = gfortran


EXE = HF


SRC = HF.f90

LIBS = -L$(LIB) -llapack -lgfortran -lrefblas

LIB = /home/borna/Desktop/Jason\ Goodpaster/Lapack/lapack-3.6.0/

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=    $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): $(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:	
	rm -f $(OBJ) $(EXE) core 
