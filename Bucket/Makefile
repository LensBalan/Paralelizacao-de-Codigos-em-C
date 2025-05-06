CC = gcc
OPENMPFLAGS = -fopenmp
MPICC = mpicc

#Executaveis
EXEC_SEQ = programa_seq
EXEC_PAR = programa_OpenMP
EXEC_MPI = programa_MPI

#Arqs fonte
SRC_SEQ = main.c
SRC_PAR = main_paralelo_OpenMP.c
SRC_MPI = main_paralelo_MPI.c


all: $(EXEC_SEQ) $(EXEC_PAR) $(EXEC_MPI)

#sequencial
$(EXEC_SEQ): $(SRC_SEQ)
	$(CC) -o $(EXEC_SEQ) $(SRC_SEQ)

#OpenMP
$(EXEC_PAR): $(SRC_PAR)
	$(CC) $(OPENMPFLAGS) -o $(EXEC_PAR) $(SRC_PAR)

#MPI
$(EXEC_MPI): $(SRC_MPI)
	$(MPICC) -o $(EXEC_MPI) $(SRC_MPI)
