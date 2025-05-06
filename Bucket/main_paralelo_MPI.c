#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX 10  // Número de elementos em cada bucket

// Função para encontrar o valor máximo em um array
float findMax(float array[], int n) {
    float max = array[0];
    for (int i = 1; i < n; i++) {
        if (array[i] > max) {
            max = array[i];
        }
    }
    return max;
}

// Função para inserir um elemento no balde
void insertionSort(float bucket[], int n) {
    for (int i = 1; i < n; i++) {
        float key = bucket[i];
        int j = i - 1;
        while (j >= 0 && bucket[j] > key) {
            bucket[j + 1] = bucket[j];
            j--;
        }
        bucket[j + 1] = key;
    }
}

// Função para realizar o bucket sort
void bucketSort(float array[], int n, int rank, int size) {
    // Encontrar o valor máximo no array
    float max = findMax(array, n);

    // Inicializar os baldes
    int bucketCount = MAX;
    float buckets[bucketCount][n];
    int bucketSizes[bucketCount];

    // Inicializar os tamanhos dos buckets
    for (int i = 0; i < bucketCount; i++) {
        bucketSizes[i] = 0;
    }

    // Distribuir os elementos nos baldes
    for (int i = 0; i < n; i++) {
        int index = (array[i] * bucketCount) / (max + 1);
        buckets[index][bucketSizes[index]++] = array[i];
    }

    // Ordenar cada balde individualmente
    for (int i = 0; i < bucketCount; i++) {
        if (bucketSizes[i] > 0) {
            insertionSort(buckets[i], bucketSizes[i]);
        }
    }

    // Concatenar os baldes ordenados de volta no array original
    int index = 0;
    for (int i = 0; i < bucketCount; i++) {
        for (int j = 0; j < bucketSizes[i]; j++) {
            array[index++] = buckets[i][j];
        }
    }
}

// Função para imprimir o array
void printArray(float array[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%.2f ", array[i]);
    }
    printf("\n");
}

// Função para ler os dados de um arquivo
int readFile(const char *filename, float **array) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        return -1;
    }

    int n = 0;
    while (fscanf(file, "%f", &(*array)[n]) == 1) {
        n++;
        *array = realloc(*array, sizeof(float) * (n + 1));
    }
    fclose(file);
    return n;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    float *array = malloc(sizeof(float) * 1);
    const char *filename = "ummilhao.txt";  // Substitua pelo nome do seu arquivo

    int n = readFile(filename, &array);
    if (n == -1) {
        MPI_Finalize();
        return 1;
    }

    //Dividi o array entre os processos
    int local_n = n / size;
    float *local_array = malloc(sizeof(float) * local_n);
    MPI_Scatter(array, local_n, MPI_FLOAT, local_array, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    
    double start_time = MPI_Wtime();

    //Ordena localmente
    bucketSort(local_array, local_n, rank, size);

    //Reuni os arrays ordenados
    float *sorted_array = NULL;
    if (rank == 0) {
        sorted_array = malloc(sizeof(float) * n);
    }
    MPI_Gather(local_array, local_n, MPI_FLOAT, sorted_array, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //Ordena o array final no processo 0, evitar sobreposicao
    if (rank == 0) {
        bucketSort(sorted_array, n, rank, size);

        double end_time = MPI_Wtime();
        printf("Array ordenado: \n");
        printArray(sorted_array, n);

        printf("Tempo de execução: %.6f segundos\n", end_time - start_time);

        free(sorted_array);
    }

    free(array);
    free(local_array);

    MPI_Finalize();
    return 0;
}
