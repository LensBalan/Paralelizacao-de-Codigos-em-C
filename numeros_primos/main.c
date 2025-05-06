#include <stdio.h>
#include <math.h>
#include <time.h> 

#define MAX 30000000

int main() {
    
    clock_t start = clock();

    /* Create an array of values, where '1' indicates that a number is prime.
     * Start by assuming all numbers are prime by setting them to 1.
     */
    char primes[MAX];
    for (long i = 0; i < MAX; i++) {
        primes[i] = 1;
    }

    /* Loop through a portion of the array (up to the square root of MAX). If
     * it's a prime, ensure all multiples of it are set to zero (false), as they
     * clearly cannot be prime.
     */
    int limit = sqrt(MAX) + 1;
    for (int i = 2; i < limit; i++) {
        if (primes[i - 1]) {
            for (long j = i * i; j <= MAX; j += i) {
                primes[j - 1] = 0;
            }
        }
    }

    /* Output the results */
    int count = 0;
    for (long i = 2; i <= MAX; i++) {
        if (primes[i - 1]) {
            // printf("%d\n", i);
            count++;
        }
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    printf("There were %d primes up to %d\n", count, MAX);
    printf("Tempo de execução: %.6f segundos\n", time_spent);

    return 0;
}