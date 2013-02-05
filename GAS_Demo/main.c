/*
	Sample demo software for Genetic Algorithms Homework
	Antenna design using GAs.
	The task is to find an optimal wire antenna design (shape) knowing the ranges of the beamwidths in azimuth and elevation
	and the target gain.

	Problem statement:
	Using GAs design a 3D wire antenna which will have a target gain of T (dB) at a certain frequency
	given the beamwidth values of the azimuth and elevation angles, a fixed number of segments and a fixed segment size.
	The output of the application should be a vector of coordinates for the segment points and the error from the desired target gain.

	Input to the application:
	  - the maximum length of the antenna in segments (segments = points-1)
	  - the range of the elevation
	  - the range of the azimuth
	  - the maximum segment size
	  - target gain
	  - frequency

	Output:
	  - coordinates of the antenna segments
	  - error between output gain and desired gain

	Hints:

	  - the maximum segment len is given by:
		Lmax = 1.5*lambda     (lambda - wavelength at the given frequency)

		lambda = C/f	      (C - lightspeed [m/s] and f - frequency [Hz = 1/s] )

	  - you have to split the segment length interval which ranges [0, Lmax]
	    so that you can get the proper binary representation

	    e.g.:
		  if Lmax = 9 at f = 1600MHz => sub_intervals = Lmax/(1.5*lambda) = 32  =  2^5 => n = 5 bits / gene

		    lambda = C/f

	  - as we are considering a 3D antenna we have 3 coordinates for each point (Px,Py,Pz)

	      e.g.:
		    for 7 segments we have 8 points each with 3 coordintes => 24 genes in each chromosome
		    but we have 5 bits per gene in the representation so we get 120 bits per chromosome - representation size

	  - as we want to get a closer gain to the target one the fitness function is given by:

	      f = sum_over_all_theta_phi(gain(theta, phi) - T)^2

	      where,

		T - target gain
		G(theta, phi) - gain for a given beamwidth of azimuth and elevation

		  G(theta, phi) = 16/(sin(theta)*sin(phi))  - ellipsoid antenna pattern of activity


	To be used as validator in the CI class @ TUM WSS 2012-13
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/* GAs parameters */
#define POPULATION_SIZE     50            // chromosomes
#define MAX_GENERATIONS     150         // number of generations to evolve
#define XOVER_PROB          0.8           // crossover probability
#define MUTATION_PROB       0.15          // mutation probability
#define MAX_INPUT_PTS	    1000

/////////////////////////////// PROBLEM SPECIFIC CODE //////////////////////////////////


/* upper and lower interval for segment size*/
double Lmin = 0.0f;
double Lmax = 0.0f;
int representation_size=0; 		  // bits
/* number of segments of the antenna */
int nr_seg = 0.0f;
/* the range of the azimuth angle */
int phi_min = 0;
int phi_max = 0;
/* the range of the elevation angle */
int theta_min = 0;
int theta_max = 0;
/* the target gain */
double target_gain = 0.0f;
/* frequency value */
int frequency = 0;  //Mhz
/* number of data points */
int data_pts = 0;

/* log data */
FILE *fp;

/* compute the gain */
double g(int phi, int theta){
  return (double)(16.0/(sin(phi)*sin(theta)));
}

/* the fitness function to be minimised - minimum defines the gain closest to the desired one */
double f(int phi_min, int phi_max, int theta_min, int theta_max){
  double ret = 0.0f;
  for(int i=phi_min;i<phi_max;++i){
      for(int j=theta_min;j<theta_max;++j){
            ret += pow(g(i, j) - target_gain, 2);
        }
    }
  return ret;
}

/////////////////////////////END OF PROBLEM SPECIFIC CODE /////////////////////////////


//////////////////////////////////// GAS SPECIFIC CODE ///////////////////////////////////////////

/* chromosome abstraction */
struct chromosome{
  /* binary chromosome representation */
  int *b;
  /* converted binary to decimal value encoded */
  double x_star;
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

/* random value generator within the problem bounds */
double randomize(double l, double h)
{
  return (((double)(rand()%1000)/1000.0)*(h - l) + l);
}

/* mapping / encoding function from binary to decimal given a binary string */
double encode_binary_chromosome(int *b)
{
  double ret = 0.0f;
  for(int i=0;i<representation_size; ++i){
      ret+=b[i]*pow(2, i);
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
  /* compute the decimal representation of the genes */
  c->x_star = encode_binary_chromosome(c->b); //randomize(Lmin, Lmax);
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
  for(int i=0; i<p->size; i++){
      init_chromosome(&p->c[i]);
    }
}

/* evaluate function, takes a user defined function and computes it for every chromosome */
void evaluate_population(struct population *p)
{
  for(int i=0; i<p->size; i++){
      p->c[i].fitness = f(p->c[i].x_star);
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
  p->c[POPULATION_SIZE].x_star = p->c[p->best_chromosome_idx].x_star;
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
      if(p->c[i].fitness > p->c[i+1].fitness){
          if(p->c[i].fitness >= best->fitness){
              best->fitness = p->c[i].fitness;
              best_idx = i;
            }
          if(p->c[i+1].fitness <= worst->fitness){
              worst->fitness = p->c[i+1].fitness;
              worst_idx = i+1;
            }
        }
      else{
          if(p->c[i].fitness <= worst->fitness){
              worst->fitness = p->c[i].fitness;
              worst_idx = i;
            }
          if(p->c[i+1].fitness >= best->fitness){
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
  if(best->fitness >= p->c[POPULATION_SIZE].fitness){
      p->c[POPULATION_SIZE].x_star = p->c[best_idx].x_star;
      p->c[POPULATION_SIZE].fitness = p->c[best_idx].fitness;
    }
  else{
      p->c[worst_idx].x_star = p->c[POPULATION_SIZE].x_star;
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
  init_chromosome(p1);
  /* cross over loop */
  for(int i=0; i< p->size; ++i){
      prob_xover = rand()%1000/1000.0;
      if(prob_xover < XOVER_PROB){candidate
          cnt++;
          if(cnt%2==0){
              double tmp;
              tmp = p1->x_star;
              p1->x_star = p->c[i].x_star;
              p->c[i].x_star = tmp;
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
          p->c[i].x_star = randomize(Lmin, Lmax);
        }
    }
}

/* print the state of the current evolution generation */
void report_state(struct population *p)
{
  printf("Generation: %d | Best fitness: %lf\n", p->gen, p->c[POPULATION_SIZE].fitness);
  fprintf(fp, "%d,%lf\n", p->gen, p->c[POPULATION_SIZE].fitness);
}

/* entry point */
int main(int argc, char* argv[]){
  fp = fopen("gas_demo_log.txt","w+");
  srand(time(NULL));
  printf("\n\nSimulation for GAs started...\n\n");
  struct population *p = (struct population*)calloc(1, sizeof(struct population));
  struct population *newp = (struct population*)calloc(1, sizeof(struct population));
  init_population(p, POPULATION_SIZE);
  init_population(newp, POPULATION_SIZE);
  evaluate_population(p);
  select_best(p);
  report_state(p);
  while(p->gen < MAX_GENERATIONS ){
      p->gen++;
      apply_selection(p, newp);
      apply_crossover(p);
      apply_mutation(p);
      report_state(p);
      evaluate_population(p);
      apply_elitism(p);
    }
  printf("\nEvolution is completed...\n\n");
  printf("\nBest chromosome: %lf | Best fitness: %lf\n\n", p->c[POPULATION_SIZE].x_star, p->c[POPULATION_SIZE].fitness);
  printf("\nSimulation ended.\n\n");
  fclose(fp);
  return EXIT_SUCCESS;
}
