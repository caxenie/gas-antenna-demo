/*
    Sample demo software for Genetic Algorithms Homework
    Antenna design using GAs.
    The task is to find an optimal wire antenna design (shape) knowing the number of points,
    the frequency and the desired gain.

    Problem statement:
    Using GAs design a 3D wire antenna which will have a target gain of T (dB) at a certain frequency, f
    given a fixed number of points or segments.
    The output of the application should be a vector of (x, y, z) coordinates for each point and the error
    from the desired target gain.

    Input to the application:
      - the maximum length of the antenna in points (segments = points-1)
      - frequency
      - target gain

    Output:
      - coordinates of the antenna points
      - error between output gain and desired gain

    Hints:

      - the maximum segment len is given by:
        Lmax = 0.5*lambda     (lambda - wavelength at the given frequency (Hz), Lmax in m)

        lambda = C/f	      (C - lightspeed [m/s] and f - frequency [Hz = 1/s] )

      - each point is corresponding to a gene

      - as we are considering a 3D antenna we have 3 coordinates for each point (Px,Py,Pz)

          e.g.:
            for 8 points we have 7 segments and 8 genes in each chromosome
            but we have 5 bits per coordinate in the representation so we get 3x5 = 15bits/gene and 120 bits per chromosome - representation size

      - as we want to get a closer gain to the target one the fitness function is given by:

          f = gain(theta, phi) - T

          where,

        T - target gain
        G(theta, phi) - gain for a given beamwidth of azimuth and elevation of the antenna

          G(theta, phi) = 16/(sin(theta)*sin(phi))  - ellipsoid antenna pattern of activity


    To be used as validator in the CI class @ TUM WSS 2012-13
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/* GAs parameters */
#define POPULATION_SIZE     10            // chromosomes
#define MAX_GENERATIONS     10000         // number of generations to evolve
#define XOVER_PROB          0.7           // crossover probability
#define MUTATION_PROB       0.25          // mutation probability
#define C_AIR               (double)300000000        // m/s
#define LAMBDA(f)           (double)C_AIR/f           // wavelength for a given frequency
#define TO_DEGREES(x)        (double)(x*(180/M_PI))
#define RAND_BOUNDED        0.123
#define MAX_FITNESS         (double)9999.9999
#define BITS_PER_GENE       15

/////////////////////////////// PROBLEM SPECIFIC CODE //////////////////////////////////


/* upper and lower interval for segment size*/
double Lmax = 0;
int representation_size=0; 		  // bits in the chromosome
int gene_size= 0;                 // bits / gene (3 coordinates: x, y, z)
/* number of points and segments of the antenna */
int nr_points = 0;
/* the target gain */
double target_gain = 0.0f;
/* frequency value */
int frequency = 0;  //Mhz
/* number of data points */
int data_pts = 0;

/* log data */
FILE *fp;

/* struct to embed the points in the wired antenna design */
struct point{
    double px;
    double py;
    double pz;
};

/* compute the azimuth and elevation angles of a segment from the coordinates */
double * compute_angles(struct point *start, struct point *end){
    double *angles = (double*)calloc(2, sizeof(double));
    double r = sqrt(pow(end->px - start->px, 2) + pow(end->py - start->py, 2) + pow(end->pz - start->pz, 2));
    double theta = acos(fabs(end->pz - start->pz)/r);
    if((int)(theta)>=1){
        theta = RAND_BOUNDED;
    }
    double phi = atan(fabs(end->py - start->py)/fabs(end->px - start->px));
    if((int)(phi)==0){
        phi = RAND_BOUNDED;
    }
    angles[0] = theta;
    angles[1] = phi;
    return angles;
}

/* compute the gain of an antenna given by the 2 angles for an elliptic antenna */
double compute_gain(double theta, double phi){
    double ret = 0.0f;
    ret = (double)(16.0/(sin(phi)*sin(theta)));
    ret = 10*log10(ret); // convert to dB
    return ret;
}

/* the fitness function to be minimised - minimum defines the gain closest to the desired one */
double f(double phi, double theta){
    double ret = 0.0f;
    ret = compute_gain(phi, theta) - target_gain;
    return ret;
}

/////////////////////////////END OF PROBLEM SPECIFIC CODE /////////////////////////////


//////////////////////////////////// GAS SPECIFIC CODE ///////////////////////////////////////////

/* chromosome abstraction */
struct chromosome{
    /* binary chromosome representation */
    int *b;
    /* converted binary to decimal point coords */
    struct point* pts;
    /* the fitness value of the chromosome */
    double fitness;
    /* relative fitness */
    double rfitness;
    /* cumulative fitness */
    double cfitness;
};

/* chromosomes population */
struct population{
    /* size in chromosomes */
    int size;
    /* chromosomes in the population */
    struct chromosome *c;
    /* current evolution generation */
    int gen;
    /* index of the fittest chromosome */
    int best_chromosome_idx;
};

double randomize(double l, double h)
{
    return (((double)(rand()%1000)/1000.0)*(h - l) + l);
}

/* mapping / decoding function from binary to decimal given a binary string */
struct point* decode_binary_chromosome(int *b)
{
    struct point* ret = (struct point*) calloc(nr_points, sizeof(struct point));
    for(int j = 0;j<nr_points;++j){
        for(int t = 0;t<5/* bits / coord */;++t){
            ret[j].px += (double)(b[t*j]*pow(2, t));
            ret[j].py += (double)(b[t*j+5]*pow(2, t));
            ret[j].pz += (double)(b[t*j+10]*pow(2, t));
        }
    }
    return ret;
}

/* initialize a chromosome */
void init_chromosome(struct chromosome *c){
    /* create a binary chromosome representation with random genes */
    c->b = (int*)calloc(representation_size, sizeof(int));
    for(int i=0;i<representation_size;++i){
        c->b[i] = ((rand()%1000/1000.0)<0.5)?0:1;
    }
    /* compute the decimal representation of the genes - get the real valued coords */
    c->pts = decode_binary_chromosome(c->b);
    /* fitness values */
    c->fitness = 0.0f;
    c->rfitness = 0.0f;
    c->cfitness = 0.0f;
}

/* initialize a chromosomes population with given parameters */
void init_population(struct population *p, int psize){
    p->size = psize;
    p->gen = 0;
    p->best_chromosome_idx = 0;
    p->c = (struct chromosome*)calloc(p->size, sizeof(struct chromosome));
    for(int i=0; i<p->size + 1  ; i++){// one more position to store the best at the end
        init_chromosome(&p->c[i]);
    }
}

/* evaluate function, takes a user defined function and computes it for every chromosome */
void evaluate_population(struct population *p)
{
    for(int i=0; i<p->size; i++){
        for(int j=1;j<nr_points;++j){
            double* angles = compute_angles(&(p->c[i].pts[j-1]), &(p->c[i].pts[j]));
            p->c[i].fitness = f(angles[0], angles[1]);
            if(isnan(p->c[i].fitness)!=0 || isinf(p->c[i].fitness)!=0)
                p->c[i].fitness = MAX_FITNESS;
        }
    }
}

/* select the best (fittest) chromosome in the population */
void select_best(struct population *p)
{
    p->best_chromosome_idx = 0;
    for(int i=0; i<p->size; ++i){
        /* the last entry in the population is the best chromosome */
        if(p->c[i].fitness > p->c[POPULATION_SIZE].fitness){
            p->best_chromosome_idx = i;
            p->c[POPULATION_SIZE].fitness = p->c[i].fitness;
        }
    }
    /* found the fittest then copy the genes */
    for(int i=0;i<nr_points;++i){
        p->c[POPULATION_SIZE].pts[i].px = p->c[p->best_chromosome_idx].pts[i].px;
        p->c[POPULATION_SIZE].pts[i].py = p->c[p->best_chromosome_idx].pts[i].py;
        p->c[POPULATION_SIZE].pts[i].pz = p->c[p->best_chromosome_idx].pts[i].pz;
    }
}

/* apply elitism so that if the previous best chromosome is better than the
 * current generation best the first will replace the worst chromosome in the
 * current generation.
 */
void apply_elitism(struct population *p)
{
    struct chromosome *best = (struct chromosome*)calloc(1, sizeof(struct chromosome));
    struct chromosome *worst= (struct chromosome*)calloc(1, sizeof(struct chromosome));
    int best_idx = 0, worst_idx = 0;
    init_chromosome(best);
    init_chromosome(worst);
    best->fitness = p->c[0].fitness;
    worst->fitness = p->c[0].fitness;

    for(int i=0;i< p->size-1;++i){
        if(p->c[i].fitness < p->c[i+1].fitness){
            if(p->c[i].fitness <= best->fitness){
                best->fitness = p->c[i].fitness;
                best_idx = i;
            }
            if(p->c[i+1].fitness >= worst->fitness){
                worst->fitness = p->c[i+1].fitness;
                worst_idx = i+1;
            }
        }
        else{
            if(p->c[i].fitness >= worst->fitness){
                worst->fitness = p->c[i].fitness;
                worst_idx = i;
            }
            if(p->c[i+1].fitness <= best->fitness){
                best->fitness = p->c[i+1].fitness;
                best_idx = i+1;
            }
        }
    }
    /* if best chromosome from the new population is better than */
    /* the best chromosome from the previous population, then    */
    /* copy the best from the new population; else replace the   */
    /* worst chromosome from the current population with the     */
    /* best one from the previous generation                     */
    if(best->fitness <= p->c[POPULATION_SIZE].fitness){
        for(int i=0;i<nr_points;++i){
            p->c[POPULATION_SIZE].pts[i].px = p->c[best_idx].pts[i].px;
            p->c[POPULATION_SIZE].pts[i].py = p->c[best_idx].pts[i].py;
            p->c[POPULATION_SIZE].pts[i].pz = p->c[best_idx].pts[i].pz;
        }
        p->c[POPULATION_SIZE].fitness = p->c[best_idx].fitness;
    }
    else{
        for(int i=0;i<nr_points;++i){
            p->c[worst_idx].pts[i].px = p->c[POPULATION_SIZE].pts[i].px;
            p->c[worst_idx].pts[i].py = p->c[POPULATION_SIZE].pts[i].py;
            p->c[worst_idx].pts[i].pz = p->c[POPULATION_SIZE].pts[i].pz;
        }
        p->c[worst_idx].fitness = p->c[POPULATION_SIZE].fitness;
    }
}

/* selection function using the elitist model in which only the
 * best chromosome survives - Winner-Take-All
 */
void apply_selection(struct population *p, struct population *newp)
{
    double sum_fit = 0.0f;
    double prob = 0.0f;
    /* find the global value of the fitness of the population */
    for(int i=0; i< p->size; ++i){
        sum_fit+=p->c[i].fitness;
    }
    /* compute the relative fitness of the population */
    for(int i=0; i<p->size; ++i){
        p->c[i].rfitness = p->c[i].fitness/sum_fit;
    }
    p->c[0].cfitness = p->c[0].rfitness;

    /* compute the cumulative fitness of the population */
    for(int i=1; i< p->size; ++i){
        p->c[i].cfitness = p->c[i-1].cfitness + p->c[i].rfitness;
    }
    /* select the survivors using the cumulative fitness */
    for(int i=0;i<p->size;++i){
        prob = rand()%1000/1000.0;
        if(prob < p->c[0].cfitness)
            newp->c[i] = p->c[0];
        else
        {
            for(int j=0; j<p->size; ++j){
                if(prob>=p->c[j].cfitness && prob<p->c[j+1].cfitness)
                    newp->c[i] = p->c[j+1];
            }
        }
    }
    /* one the new population is created copy it back in the working var */
    for(int i=0 ;i<p->size; ++i)
        p->c[i] = newp->c[i];
}

/* apply the single point crossover operator which takes 2 parents */
void apply_crossover(struct population *p)
{
    /* counter of members chosen */
    int cnt = 0;
    /* probability to xover */
    double prob_xover = 0.0f;
    /* the two parent containers init */
    struct chromosome *p1 = (struct chromosome*)calloc(1, sizeof(struct chromosome));
    struct point *tmp = (struct point*)calloc(1, sizeof(struct point));
    init_chromosome(p1);
    /* cross over loop */
    for(int i=0; i< p->size; ++i){
        prob_xover = rand()%1000/1000.0;
        if(prob_xover < XOVER_PROB){
            cnt++;
            if(cnt%2==0){
                for(int j=0;j<nr_points;++j){
                    // swap
                    tmp->px = p1->pts[j].px;
                    tmp->py = p1->pts[j].py;
                    tmp->pz = p1->pts[j].pz;

                    p1->pts[j].px = p->c[i].pts[j].px;
                    p1->pts[j].py = p->c[i].pts[j].py;
                    p1->pts[j].pz = p->c[i].pts[j].pz;

                    p->c[i].pts[j].px = tmp->px;
                    p->c[i].pts[j].py = tmp->py;
                    p->c[i].pts[j].pz = tmp->pz;
                }
            }
            else
            {
                p1 = &p->c[i];
            }
        }
    }
}

/* apply mutation - random uniform mutation of the genes */
void apply_mutation(struct population *p)
{
    double prb = 0.0f;
    for(int i=0;i<p->size;++i){
        prb = rand()%1000/1000.0;
        if(prb < MUTATION_PROB){
            for(int j=0;j<nr_points;++j){
                p->c[i].pts[j].px = randomize(0, Lmax);
                p->c[i].pts[j].py = randomize(0, Lmax);
                p->c[i].pts[j].pz = randomize(0, Lmax);
            }
        }
    }
}

/* print the state of the current evolution generation */
void report_state(struct population *p)
{
    printf("Generation: %d | Best fitness: %lf\n", p->gen, p->c[POPULATION_SIZE].fitness);
    for(int i=0;i<nr_points;++i){
        fprintf(fp, "%d,%lf\n", p->gen, p->c[POPULATION_SIZE].fitness);
    }
}

/* get the input data and store it locally for the population initialization */
void get_input_data(){
    /* loop and get training data */
    while(scanf("%d,%d,%lf\n", &nr_points, &frequency, &target_gain)>0){
        printf("got input : Np = %d, Freq = %d, TGain = %lf \n\n", nr_points, frequency, target_gain);
    }
}

/* entry point */
int main(int argc, char* argv[]){
    fp = fopen("gas_demo_log.txt","w+");
    srand(time(NULL));
    printf("\n\nSimulation for Antenna Design GAs started...\n\n");
    struct population *p = (struct population*)calloc(1, sizeof(struct population));
    struct population *newp = (struct population*)calloc(1, sizeof(struct population));

    /* get data from the input file */
    get_input_data();
    Lmax = 0.5 * LAMBDA(frequency);                           // in m - max length of the segment
    representation_size= nr_points*BITS_PER_GENE;             // bits / chromosome
    gene_size= BITS_PER_GENE;                                 // bits / gene (3 coordinates: x, y, z)

    /* initialize the GA and start the evolution */
    init_population(p, POPULATION_SIZE);
    init_population(newp, POPULATION_SIZE);
    evaluate_population(p);
    select_best(p);
    report_state(p);
    while(p->gen < MAX_GENERATIONS){
        p->gen++;
        apply_selection(p, newp);
        apply_crossover(p);
        apply_mutation(p);
        report_state(p);
        evaluate_population(p);
        apply_elitism(p);
    }
    printf("\nEvolution is completed...\n\n");
    printf("Max segment length is %lf m\n", Lmax);
    printf("\nBest fitness: %lf\n\n", p->c[POPULATION_SIZE].fitness);
    printf("\n Best chromosome: \n");
    for(int i=0;i<nr_points;++i){
        printf(" %lf %lf %lf | ",  p->c[POPULATION_SIZE].pts[i].px, p->c[POPULATION_SIZE].pts[i].py, p->c[POPULATION_SIZE].pts[i].pz);
    }
    printf("\n");
    printf("\nSimulation ended.\n\n");
    fclose(fp);
    return EXIT_SUCCESS;
}
