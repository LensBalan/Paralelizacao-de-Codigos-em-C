#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 500000  //N de elementos a serem gerados
#define FILENAME "dados1.txt"  

int main() {
    
    srand(time(NULL));

    
    FILE *file = fopen(FILENAME, "w");
    if (file == NULL) {
        printf("Erro ao abrir o arquivo para escrita!\n");
        return 1;
    }

    
    for (int i = 0; i < N; i++) {
        float rand_num = (float)rand() / (float)RAND_MAX;  
        fprintf(file, "%.2f\n", rand_num); 
    }

 
    fclose(file);

    printf("Arquivo '%s' gerado com sucesso!\n", FILENAME);
    return 0;
}
