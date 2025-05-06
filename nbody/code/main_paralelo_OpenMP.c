#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */

#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;

long RandomSeed() {
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
    const long Q = MODULUS / MULTIPLIER;
    const long R = MODULUS % MULTIPLIER;
    long t;

    #pragma omp critical
    {
        t = MULTIPLIER * (seed % Q) - R * (seed / Q);
        if (t > 0) 
            seed = t;
        else 
            seed = t + MODULUS;
    }
    return seed;
}

double Random() {
    return ((double) RandomSeed() / MODULUS);
}

/*
 * End of the pRNG algorithm
 */

typedef struct {
    double x, y, z;
    double mass;
} Particle;

typedef struct {
    double xold, yold, zold;
    double fx, fy, fz;
} ParticleV;

void InitParticles(Particle[], ParticleV[], int);
double ComputeForces(Particle[], Particle[], ParticleV[], int);
double ComputeNewPos(Particle[], ParticleV[], int, double);

int main() {

    //N de threads a serem utilizadas
    omp_set_num_threads(12);

    double time_start, time_end; 
    Particle *particles; /* Particles */
    ParticleV *pv; /* Particle velocity */
    int npart, cnt;  /* number of times in loop */
    double sim_t;   /* Simulation time */

    fscanf(stdin, "%d %d", &npart, &cnt);

    printf("\nExecutando...\n\n");

    time_start = omp_get_wtime();
    
    particles = (Particle *) malloc(sizeof(Particle) * npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV) * npart);
    
    InitParticles(particles, pv, npart);
    sim_t = 0.0;

    while (cnt--) {
        double max_f;
        max_f = ComputeForces(particles, particles, pv, npart);
        sim_t += ComputeNewPos(particles, pv, npart, max_f);
    }

    time_end = omp_get_wtime();

    for (int i = 0; i < npart; i++)
        fprintf(stdout, "%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z);

    printf("Tempo de execução: %.6lf segundos\n", time_end - time_start);

    free(particles);
    free(pv);

    return 0;
}

void InitParticles(Particle particles[], ParticleV pv[], int npart) {
    #pragma omp parallel for
    for (int i = 0; i < npart; i++) {
        particles[i].x = Random();
        particles[i].y = Random();
        particles[i].z = Random();
        particles[i].mass = 1.0;
        pv[i].xold = particles[i].x;
        pv[i].yold = particles[i].y;
        pv[i].zold = particles[i].z;
        pv[i].fx = 0;
        pv[i].fy = 0;
        pv[i].fz = 0;
    }
}

double ComputeForces(Particle myparticles[], Particle others[], ParticleV pv[], int npart) {
    double max_f = 0.0;
    
    #pragma omp parallel for reduction(max:max_f)
    for (int i = 0; i < npart; i++) {
        double xi = myparticles[i].x, yi = myparticles[i].y;
        double fx = 0.0, fy = 0.0, rmin = 100.0;
        
        for (int j = 0; j < npart; j++) {
            double rx = xi - others[j].x, ry = yi - others[j].y;
            double mj = others[j].mass;
            double r = rx * rx + ry * ry;
            if (r == 0.0) continue;
            if (r < rmin) rmin = r;
            r = r * sqrt(r);
            fx -= mj * rx / r;
            fy -= mj * ry / r;
        }
        
        #pragma omp critical
        {
            pv[i].fx += fx;
            pv[i].fy += fy;
        }
        
        double fmag = sqrt(fx * fx + fy * fy) / rmin;
        if (fmag > max_f) {
            max_f = fmag;
        }
    }
    return max_f;
}

double ComputeNewPos(Particle particles[], ParticleV pv[], int npart, double max_f) {
    static double dt_old = 0.001, dt = 0.001;
    double dt_new;
    double a0 = 2.0 / (dt * (dt + dt_old));
    double a2 = 2.0 / (dt_old * (dt + dt_old));
    double a1 = -(a0 + a2);

    #pragma omp parallel for
    for (int i = 0; i < npart; i++) {
        double xi = particles[i].x, yi = particles[i].y;
        particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
        particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
        pv[i].xold = xi;
        pv[i].yold = yi;
        pv[i].fx = 0;
        pv[i].fy = 0;
    }

    dt_new = 1.0 / sqrt(max_f);
    if (dt_new < 1.0e-6) dt_new = 1.0e-6;
    if (dt_new < dt) {
        dt_old = dt;
        dt = dt_new;
    } else if (dt_new > 4.0 * dt) {
        dt_old = dt;
        dt *= 2.0;
    }
    return dt_old;
}
