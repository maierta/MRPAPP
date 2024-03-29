#
SHELL    = /bin/bash
SRC_DIR = ./
INCLUDES =  -I$(SRC_DIR) \
	-I$(SRC_DIR)/PartialPsimag \
	-I$(SRC_DIR)/PsimagLite \
	-I$(SRC_DIR)/PsimagLite/src \

EXENAME		     = mrpapp
CC                   = mpicxx
cc                   = mpicxx
MatrixBoundsChecking = -DNDEBUG

# FLAGS        = $(MatrixBoundsChecking)  -g  -Wall  -Werror -DUSE_MPI  -DUSE_SCGAP3D
# FLAGS        = $(MatrixBoundsChecking)  -O2 -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D
# FLAGS        = $(MatrixBoundsChecking)  -O2 -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D -DUSE_BAFEAS

FLAGS        = $(MatrixBoundsChecking)  -O2 -Warray-bounds -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D -DUSE_TBFILE
# FLAGS        = $(MatrixBoundsChecking)  -O3 -Warray-bounds -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D -DUSE_TBFILE
#FLAGS        = $(MatrixBoundsChecking)  -O3 -Warray-bounds -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D -DUSE_1BANDABWSPIN

# FLAGS        = $(MatrixBoundsChecking)  -O3 -Warray-bounds -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D -DUSE_1BAND
#FLAGS        = $(MatrixBoundsChecking)  -O2 -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D -DUSE_BILAYER_FESC
#FLAGS        = $(MatrixBoundsChecking)  -O2 -Wall  -Werror -DUSE_MPI -DUSE_SCGAP3D -DUSE_KFE2SE2

# FLAGS        = $(MatrixBoundsChecking)  -O2 -Wall  -Werror -DUSE_MPI  -DUSE_BILAYER
# FLAGS        = $(MatrixBoundsChecking)  -O2 -Wall  -Werror -DUSE_MPI

#Note that mjson does not compile with -Werror -Wall, this is

# a temporary solution :
FLAGSC	     = $(MatrixBoundsChecking)

OBJECTS      =  main.o


all: clean showEnv $(OBJECTS) $(EXENAME)

showEnv:
	@echo ---------------------------------------------------------------- Env Info:
	@echo ------------------------------------------------ Modules in use, if any
	@echo No Modules on this system.
	@echo ------------------------------------------------ Compiler Version
	$(CC) --version
	@echo ----------------------------------------------------------------
	@echo

main.o: $(SRC_DIR)/main.cpp
	$(CC) $(INCLUDES) $(FLAGS) -c $(SRC_DIR)/main.cpp

main.s: $(SRC_DIR)/main.cpp
	$(CC) -S $(INCLUDES) $(FLAGS) -c $(SRC_DIR)/main.cpp

$(EXENAME): $(OBJECTS)
	$(CC) $(INCLUDES) $(FLAGS) $(OBJECTS) -o $(EXENAME)  -lblas -llapack
	cp $(EXENAME) 	$(EXENAME)_stripped
	strip $(EXENAME)_stripped

clean:
	rm -f $(EXENAME) *.o
########################End of Makefile #################################
