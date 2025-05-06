#include <stdio.h>
#include <math.h>
#include <omp.h>

#define MAX 300000000

int main() {

    //N de threads a serem utilizadas
    omp_set_num_threads(12);

    double start_time, end_time; 

    start_time = omp_get_wtime();

    char primes[MAX];
    int limit = sqrt(MAX) + 1;

    //Inicializa o array de primos
    #pragma omp parallel for
    for (long i = 0; i < MAX; i++) {
        primes[i] = 1;
    }

    for (int i = 2; i < limit; i++) {
        if (primes[i - 1]) {
            #pragma omp parallel for
            for (long j = i * i; j <= MAX; j += i) {
                primes[j - 1] = 0;
            }
        }
    }

    //Conta e printa o N de primos encontrados
    int count = 0;
    #pragma omp parallel for reduction(+:count)
    for (long i = 2; i <= MAX; i++) {
        if (primes[i - 1]) {
            count++;
        }
    }

    end_time = omp_get_wtime();

    printf("There were %d primes up to %d\n", count, MAX);
    printf("Tempo de execução: %.6f segundos\n", end_time - start_time);

    return 0;
}