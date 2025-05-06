#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;

double Random(void) {
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
    const long Q = MODULUS / MULTIPLIER;
    const long R = MODULUS % MULTIPLIER;
    long t;

    t = MULTIPLIER * (seed % Q) - R * (seed / Q);
    if (t > 0) 
        seed = t;
    else 
        seed = t + MODULUS;
    return ((double) seed / MODULUS);
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
double ComputeForces(Particle[], Particle[], ParticleV[], int, int, int);
double ComputeNewPos(Particle[], ParticleV[], int, double);

int main(int argc, char *argv[]) {
    int npart, cnt; /* number of times in loop */
    Particle *particles; /* Particles */
    ParticleV *pv; /* Particle velocity */
    double sim_t = 0.0;  /* Simulation time */
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        scanf("%d", &npart);
        scanf("%d", &cnt);
        printf("\nExecutando...\n\n");
    }

    double start_time = MPI_Wtime();
    
    MPI_Bcast(&npart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cnt, 1, MPI_INT, 0, MPI_COMM_WORLD);


     /* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle) * npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV) * npart);

    if (rank == 0) {
        /* Generate the initial values */
        InitParticles(particles, pv, npart);
    }

    MPI_Bcast(particles, npart * sizeof(Particle), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(pv, npart * sizeof(ParticleV), MPI_BYTE, 0, MPI_COMM_WORLD);

    while (cnt--) {
        double max_f;
        max_f = ComputeForces(particles, particles, pv, npart, rank, size);
        
        MPI_Allreduce(MPI_IN_PLACE, &max_f, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        sim_t += ComputeNewPos(particles, pv, npart, max_f);
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, npart / size * sizeof(Particle), MPI_BYTE, particles, npart / size * sizeof(Particle), MPI_BYTE, MPI_COMM_WORLD);
    }

    double end_time = MPI_Wtime();

    if (rank == 0) {
        for (int i = 0; i < npart; i++) {
            printf("%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z);
        }
    }

    if (rank == 0) {
        printf("Tempo de execução: %f segundos\n", end_time - start_time);
    }

    free(particles);
    free(pv);
    MPI_Finalize();
    return 0;
}

void InitParticles(Particle particles[], ParticleV pv[], int npart) {
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

double ComputeForces(Particle myparticles[], Particle others[], ParticleV pv[], int npart, int rank, int size) {
    double max_f = 0.0;
    int i;
    int start = (rank * npart) / size;
    int end = ((rank + 1) * npart) / size;

    for (i = start; i < end; i++) {
        int j;
        double xi, yi, fx, fy, rmin = 100.0;
        xi = myparticles[i].x;
        yi = myparticles[i].y;
        fx = fy = 0.0;

        for (j = 0; j < npart; j++) {
            double rx = xi - others[j].x;
            double ry = yi - others[j].y;
            double mj = others[j].mass;
            double r = rx * rx + ry * ry;
            if (r == 0.0) continue;
            if (r < rmin) rmin = r;
            r = r * sqrt(r);
            fx -= mj * rx / r;
            fy -= mj * ry / r;
        }
        pv[i].fx += fx;
        pv[i].fy += fy;
        double force_magnitude = sqrt(fx * fx + fy * fy) / rmin;
        if (force_magnitude > max_f) max_f = force_magnitude;
    }
    
    return max_f;
}

double ComputeNewPos(Particle particles[], ParticleV pv[], int npart, double max_f) {
    int i;
    static double dt_old = 0.001, dt = 0.001;
    double a0 = 2.0 / (dt * (dt + dt_old));
    double a1 = -(a0 + 2.0 / (dt_old * (dt + dt_old)));

    for (i = 0; i < npart; i++) {
        double xi = particles[i].x;
        double yi = particles[i].y;
        particles[i].x = (pv[i].fx - a1 * xi - 2.0 / (dt_old * (dt + dt_old)) * pv[i].xold) / a0;
        particles[i].y = (pv[i].fy - a1 * yi - 2.0 / (dt_old * (dt + dt_old)) * pv[i].yold) / a0;
        pv[i].xold = xi;
        pv[i].yold = yi;
        pv[i].fx = 0;
        pv[i].fy = 0;
    }

    double dt_new = 1.0 / sqrt(max_f);
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
