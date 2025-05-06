CC = gcc
OPENMPFLAGS = -fopenmp
MPICC = mpicc
LIBS = -lm  

#executaveis
EXEC_SEQ = programa_seq
EXEC_OPENMP = programa_OpenMP
EXEC_MPI = programa_MPI

# Arqs fonte
SRC_SEQ = main.c
SRC_PAR = main_paralelo_OpenMP.c
SRC_MPI = main_paralelo_MPI.c

all: $(EXEC_SEQ) $(EXEC_OPENMP) $(EXEC_MPI)

#sequencial
$(EXEC_SEQ): $(SRC_SEQ)
	$(CC) -o $(EXEC_SEQ) $(SRC_SEQ) $(LIBS)

#OpenMP
$(EXEC_OPENMP): $(SRC_PAR)
	$(CC) $(OPENMPFLAGS) -o $(EXEC_OPENMP) $(SRC_PAR) $(LIBS)

#MPI
$(EXEC_MPI): $(SRC_MPI)
	$(MPICC) -o $(EXEC_MPI) $(SRC_MPI) $(LIBS)
