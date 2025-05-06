#!/bin/bash

#Definir a pilha como ilimitada
ulimit -s unlimited

mpirun -np 2 programa_MPI -nolocal -machinefile nodos
