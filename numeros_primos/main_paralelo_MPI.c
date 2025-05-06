#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX 300000000

int main(int argc, char **argv) {
    int rank, size;
    double start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    
    start_time = MPI_Wtime();

    //Divido o array entre os processos
    int chunk_size = MAX / size;
    int start = rank * chunk_size;
    int end = (rank == size - 1) ? MAX : start + chunk_size;

    char *local_primes = (char *)malloc(chunk_size * sizeof(char));
    if (local_primes == NULL) {
        perror("malloc");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //Init array local
    for (int i = 0; i < chunk_size; i++) {
        local_primes[i] = 1;
    }

    //Encontra os primos ate a raiz quadrada de MAX
    int limit = sqrt(MAX) + 1;

    //Array global para os primos
    char *global_primes = (char *)malloc(limit * sizeof(char));
    if (global_primes == NULL) {
        perror("malloc");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
        for (int i = 0; i < limit; i++) {
            global_primes[i] = 1;
        }

        for (int i = 2; i < limit; i++) {
            if (global_primes[i - 1]) {
                for (int j = i * i; j < limit; j += i) {
                    global_primes[j - 1] = 0;
                }
            }
        }
    }

    //Compartilha os primos encontrados com todos os processos
    MPI_Bcast(global_primes, limit, MPI_CHAR, 0, MPI_COMM_WORLD);

    for (int i = 2; i < limit; i++) {
        if (global_primes[i - 1]) {
    
            int first_multiple = (start / i) * i;
            if (first_multiple < start) {
                first_multiple += i;
            }
            if (first_multiple < i * i) {
                first_multiple = i * i;
            }

            for (int j = first_multiple; j < end; j += i) {
                if (j >= start) {
                    local_primes[j - start] = 0;
                }
            }
        }
    }

    //Sincronizacao
    MPI_Barrier(MPI_COMM_WORLD);

    //so um coleta os resultados
    if (rank == 0) {
        char *all_primes = (char *)malloc(MAX * sizeof(char));
        if (all_primes == NULL) {
            perror("malloc");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < chunk_size; i++) {
            all_primes[start + i] = local_primes[i];
        }

        //Recebe as partes dos outros processos
        for (int p = 1; p < size; p++) {
            int p_start = p * chunk_size;
            int p_end = (p == size - 1) ? MAX : p_start + chunk_size;
            MPI_Recv(&all_primes[p_start], p_end - p_start, MPI_CHAR, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        int count = 0;
        for (int i = 2; i <= MAX; i++) {
            if (all_primes[i - 1]) {
                count++;
            }
        }

        end_time = MPI_Wtime();

        printf("There were %d primes up to %d\n", count-1, MAX);
        printf("Tempo de execução: %.6f segundos\n", end_time - start_time);

        free(all_primes);
    } else {
        
        MPI_Send(local_primes, chunk_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    free(local_primes);
    free(global_primes);
    MPI_Finalize();
    return 0;
}