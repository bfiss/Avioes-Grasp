#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include "include\glpk.h"

#define AUXSIZE 20
#define MAXSIZE 55

#define SEED 0xDEAD

FILE * outFile;
/* GRASP parameters */
int alpha, maxAlpha, maxIter, maxTime, randomSeed;

/* Number of planes */
int n;

/* Plane data */
typedef struct {
	int earliest;
	int ideal;
	int latest;
	double costE;
	double costL;
	int sep[MAXSIZE];
	int pos;
} Plane;

/* Auxiliary structure to sort planes by their ideal times */
struct planeOrder {
	int ideal;
	int pos;
};

Plane planes[MAXSIZE];
struct planeOrder pOrd[MAXSIZE];

void leave(char * name);

void readInput(FILE * in);

void printResult(glp_prob * Prob, FILE * out);

/* Creates a solution in a greedy randomized way */
int createSolution( int solution[], int time[], int pos );

/* Maps the order of planes in solution to a list of constraints &ij */
void mapSolution(glp_prob * Prob, int solution[]);

/* Restriction bci = ai - bi + xi = Ti and z += ai*Ei + bi*Li */
void addBasicRestriction(glp_prob * Prob, int plane);

/* Restriction xj - xi >= Sij */
void addSeparationConstraint(glp_prob * Prob, int plane1, int plane2);

/* Restrictions xj - xi >= Sij&ij + (Li - Ej)&ji and
                xi - xj >= Sji&ji + (Lj - Ei)&ij and
				&ij + &ji = 1*/
void addOrderConstraint(glp_prob * Prob, int plane1, int plane2);

/* Swap the arrival of two adjacent planes */
void swapConstraint(glp_prob * Prob, int i, int solution[], int back);

int compIdealT ( const void *, const void * );

int main(int argc, char * argv[]) {
	int i,j;
	
	srand(SEED);

	readInput(stdin);
	
	/* Default values */
    outFile = stdout;
	maxAlpha = 2;
	maxIter = 100;
	maxTime = 30;
	randomSeed = SEED;
	/* Read arguments */
	if( argc > 6 )
		argc = 6;
	switch(argc) {
	case 6:
		if( !(randomSeed = atoi(argv[5])) )
			leave(argv[0]);
	case 5:
		if( !(maxTime = atoi(argv[4])) )
			leave(argv[0]);
	case 4:
		if( !(maxIter = atoi(argv[3])) )
			leave(argv[0]);
	case 3:
		if( !(maxAlpha = atoi(argv[2])) )
			leave(argv[0]);
	case 2:
		if( !(outFile = fopen(argv[1],"w")) )
			leave(argv[0]);
	}
	
	/* Initiate positions */
	for( i = 0 ; i < n ; ++i ) {
   		pOrd[i].ideal = planes[i].ideal;
  		pOrd[i].pos = i;
	}
	qsort (pOrd, n, sizeof(struct planeOrder), compIdealT);
	for( i = 0 ; i < n ; ++i ) {
  		planes[pOrd[i].pos].pos = i;
	}

	/* Create lp instance */
	glp_prob * Prob;
	Prob = glp_create_prob();
	glp_set_prob_name(Prob, "Airplane Landing Problem");
	glp_set_obj_name(Prob, "Cost");
	
	/* Create basic constraints */
	for( i = 0 ; i < n ; ++i ) {
        addBasicRestriction(Prob,i);
	}
	
	glp_create_index(Prob);
	
	/* Create separation constraints and order variables (&ij) if necessary */
	for( i = 0 ; i < n ; ++i ) {
		for( j = i+1 ; j < n ; ++j ) {
			if( planes[i].latest >= planes[j].earliest &&
			    planes[j].latest >= planes[i].earliest ) {
                addOrderConstraint(Prob,i,j);
			} else if ( planes[i].latest < planes[j].earliest &&
						planes[i].latest + planes[i].sep[j] >= planes[j].earliest ) {
                addSeparationConstraint(Prob, i, j);
			} else if ( planes[j].latest < planes[i].earliest &&
						planes[j].latest + planes[j].sep[i] >= planes[i].earliest ) {
                addSeparationConstraint(Prob, j, i);
			}
		}
	}

	/* Write problem in MPS format so glpsol can (try to) solve it */
	glp_write_mps(Prob, GLP_MPS_FILE, NULL,"mpsProblem.txt");
	
	glp_delete_index(Prob);
	glp_create_index(Prob);
	
	/* GRASP */
	
	/* Data to handle glp solving, time checking and solution generating */
	glp_smcp * param = malloc(sizeof(glp_smcp));
	glp_init_smcp(param);
	param->msg_lev = GLP_MSG_ERR;
	int solution[MAXSIZE], timeAux[MAXSIZE], t;
	double currResult = DBL_MAX, bestResult = DBL_MAX;
	alpha = 0;
	time_t start, curr;
	time(&start);
	
	for( t = 0 ; t < maxIter ; ++t ) {
		/* Greedy solution generation */
		while(createSolution(solution,timeAux,0))
			alpha = n;
		
		/* Building the right constraints */
		mapSolution(Prob,solution);
		
		/* Solving with glpsol */
		param->presolve = GLP_ON;
		glp_simplex(Prob,param);
		param->presolve = GLP_OFF;
		currResult = glp_get_obj_val(Prob);
		
		/* Local search using the first increase */
		for( i = 0 ; i < n-1 ; ++i ) {

			/* Swap two adjacent planes */
			swapConstraint(Prob,i,solution,0);
			glp_simplex(Prob,param);
			
			/* Check for improvements */
			if( GLP_OPT == glp_get_status(Prob) && glp_get_obj_val(Prob) < currResult ) {
				
				currResult = glp_get_obj_val(Prob);
				
				/* Changing the solution */
				int swp;
				swp = solution[i];
				solution[i] = solution[i+1];
				solution[i+1] = swp;
				
				/* Restarting */
				i = -1;
			} else
				swapConstraint(Prob,i,solution,1);
		}
		
		/* Checking improvements */
		if( bestResult > currResult ) {
		    bestResult = currResult;
		    for( i = 0 ; i < n ; ++i )
				planes[solution[i]].pos = i;
		}
		
		/* Choosing alpha */
		alpha = rand()%(maxAlpha+1);
		
		/* Is our time up? */
		time(&curr);
		if( difftime(curr,start) > maxTime )
		    break;
	}
	
	/* Print Answer */
	printResult(Prob, stdout);
	if( outFile ) {
		printResult(Prob, outFile);
		fclose(outFile);
	}

	return 0;
}

void leave(char * name) {
	printf("Usage: %s [output_file [maximum_alpha [iteration_limit [time_limit [random_seed]]]]]\n",name);
	exit(1);
}

void readInput(FILE * in) {
	int i,j;
	
	fscanf(in,"%i",&n);
	fscanf(in,"%i",&i);
	for( i = 0 ; i < n ; ++i ) {
		fscanf(in,"%i",&j);
		fscanf(in,"%i",&(planes[i].earliest));
		if( j > planes[i].earliest )
		    planes[i].earliest = j;
		fscanf(in,"%i",&(planes[i].ideal));
		fscanf(in,"%i",&(planes[i].latest));
		fscanf(in,"%lf",&(planes[i].costE));
		fscanf(in,"%lf",&(planes[i].costL));
		for( j = 0 ; j < n ; ++j )
	        fscanf(in,"%i",&(planes[i].sep[j]));
	}
}

void printResult(glp_prob * Prob, FILE * out) {
	int i;
	char buf[AUXSIZE];
	glp_smcp * param = malloc(sizeof(glp_smcp));
	glp_init_smcp(param);
	param->msg_lev = GLP_MSG_ERR;
	param->presolve = GLP_ON;
	int solution[MAXSIZE];

	for( i = 0 ; i < n ; ++i )
		solution[planes[i].pos] = i;

	mapSolution(Prob,solution);
	glp_simplex(Prob,param);
	
	fprintf(out,"Best found solution's value: %lf\n\n",glp_get_obj_val(Prob));
	
	double time;
	for( i = 0 ; i < n ; ++i ) {
        sprintf(buf,"x%i",solution[i]);
		time = glp_get_col_prim(Prob, glp_find_col(Prob, buf));
		fprintf(out,"The %i-th airplane to arrive is airplane %i, at the time %lf\n",
		        i+1,solution[i]+1,time);
	}
}

/* Creates a solution in a greedy randomized way */
int createSolution( int solution[], int time[], int pos ) {
	int i,j;
    int valid[MAXSIZE];
    double cost[MAXSIZE];

    /* Initialize valid[] */
    for( i = 0 ; i < n ; ++i )
		valid[i] = 1;

	/* Find out who has been included */
	for( i = 0 ; i < pos ; ++i )
	    valid[solution[i]] = 0;

    /* For all planes not yet included, find the one
	   whose latest arrival is minimal. */
	int arrivalLimit = INT_MAX;
	for( i = 0 ; i < n ; ++i )
		if( valid[i] && arrivalLimit > planes[i].latest )
			arrivalLimit = planes[i].latest;

	/* Find all the planes that can be placed next, their
	   minimum arrival times and respective costs */
	int arrival = 0;
	double posDiff;
	for( i = 0 ; i < n ; ++i )
		if( valid[i] ) {
			/* Respect the time distance between planes */
	    	for( j = pos - 1 ; j >= 0 ; --j )
				if( arrival < planes[solution[j]].sep[i] + time[solution[j]] )
			    	arrival = planes[solution[j]].sep[i] + time[solution[j]];
			/* Respect own plane limits */
			if( arrival > planes[i].latest )
		    	return 1; /* Impossible */
			if( arrival < planes[i].earliest )
		    	arrival = planes[i].earliest;
			/* Check if this plane has to arrive after another one */
			if( arrival > arrivalLimit ) {
				valid[i] = 0;
				continue;
			}
			/* Calculate time and cost of arriving the plane now */
			time[i] = arrival;
			posDiff = (pos - planes[i].pos);
			if( posDiff < 0 )
			    posDiff = -posDiff;
			cost[i] = posDiff;
		}

	/* Find out the best option */
	double bestOption;
tryAgain:
	bestOption = DBL_MAX - 2* alpha;

	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && bestOption > cost[i] )
			bestOption = cost[i];

	/* Find out the number of options */
	int numValid;
	numValid = 0;
	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && cost[i] <= bestOption + alpha )
			++numValid;

	if( !numValid )
	    return 1;

	/* Select an option close or equal to the best */
	int option, counter;
	option = rand()%numValid + 1;
	counter = 0;
	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && cost[i] <= bestOption + alpha )
	        if( ++counter == option )
	            break;
	solution[pos] = i;

	/* Check if the solution is done */
	if( ++pos == n )
	    return 0;

	/* Check if the rest of the solution is feasible */
    if( createSolution( solution, time, pos ) ) {
		valid[i] = 0;
		--pos;
		goto tryAgain;
	}

	return 0;
}

/* Maps the order of planes in solution to a list of constraints &ij */
void mapSolution(glp_prob * Prob, int solution[]) {
	int i,j;
	char buf[AUXSIZE];
	
	int uij;
	for( i = 0 ; i < n ; ++i )
		for( j = i+1 ; j < n ; ++j ) {
			sprintf(buf,"u%i,%i",solution[i],solution[j]);
			if( (uij = glp_find_col(Prob, buf)) ) {
                glp_set_col_bnds(Prob, uij, GLP_FX, 1, 1);
                glp_set_col_kind(Prob, uij, GLP_CV);
                sprintf(buf,"u%i,%i",solution[j],solution[i]);
				glp_set_col_bnds(Prob, uij=glp_find_col(Prob, buf), GLP_FX, 0, 0);
                glp_set_col_kind(Prob, uij, GLP_CV);
			}
		}
}

/* Restriction bci = ai - bi + xi = Ti and z += ai*Ei + bi*Li */
void addBasicRestriction(glp_prob * Prob, int plane) {
    int cardinal, constr[4], i = plane, cardRow;
	double cValues[4];
	char buf[AUXSIZE];

	cardinal = glp_add_cols(Prob, 3);
	
	sprintf(buf,"a%i",i);
	glp_set_col_name(Prob, cardinal, buf);
	glp_set_col_bnds(Prob, cardinal, GLP_LO, 0, 0);
	glp_set_obj_coef(Prob, cardinal, planes[i].costE);
	
	sprintf(buf,"b%i",i);
	glp_set_col_name(Prob, cardinal+1, buf);
	glp_set_col_bnds(Prob, cardinal+1, GLP_LO, 0, 0);
	glp_set_obj_coef(Prob, cardinal+1, planes[i].costL);
	
	sprintf(buf,"x%i",i);
	glp_set_col_name(Prob, cardinal+2, buf);
	if( planes[i].earliest == planes[i].latest )
	    glp_set_col_bnds(Prob, cardinal+2, GLP_FX, planes[i].earliest, 0);
	else
		glp_set_col_bnds(Prob, cardinal+2, GLP_DB, planes[i].earliest, planes[i].latest);
	
    cardRow = glp_add_rows(Prob, 1);
    
	sprintf(buf,"bc%i",i);
	glp_set_row_name(Prob, cardRow, buf);
	glp_set_row_bnds(Prob, cardRow, GLP_FX, planes[i].ideal, 0);
	
	constr[3] = 1 + (constr[2] = 1 + (constr[1] = cardinal));
	cValues[3] = cValues[1] = 1;
	cValues[2] = -1;
	
	glp_set_mat_row(Prob, cardRow, 3, constr, cValues);
}

/* Restriction xj - xi >= Sij */
void addSeparationConstraint(glp_prob * Prob, int plane1, int plane2) {
    int cardinal, constr[3], i = plane1, j = plane2;
	double cValues[3];
	char buf[AUXSIZE];
	
    cardinal = glp_add_rows(Prob, 1);
    
	sprintf(buf,"S%i,%i",i,j);
	glp_set_row_name(Prob, cardinal, buf);
	glp_set_row_bnds(Prob, cardinal, GLP_LO, planes[i].sep[j], 0);
	
	sprintf(buf,"x%i",j);
	constr[1] = glp_find_col(Prob, buf);
	
	sprintf(buf,"x%i",i);
	constr[2] = glp_find_col(Prob, buf);
	
	cValues[1] = 1;
	cValues[2] = -1;
	
	glp_set_mat_row(Prob, cardinal, 2, constr, cValues);
}

/* Restrictions xj - xi >= Sij&ij + (Li - Ej)&ji and
                xi - xj >= Sji&ji + (Lj - Ei)&ij and
				&ij + &ji = 1*/
void addOrderConstraint(glp_prob * Prob, int plane1, int plane2) {
    int cardinal, constr[3], i = plane1, j = plane2;
	double cValues[3];
	char buf[AUXSIZE];
	int xi, xj, uij;
	
	sprintf(buf,"x%i",i);
	xi = glp_find_col(Prob, buf);
    sprintf(buf,"x%i",j);
	xj = glp_find_col(Prob, buf);
	
	uij = glp_add_cols(Prob, 2);
	
	sprintf(buf,"u%i,%i",i,j);
	glp_set_col_name(Prob, uij, buf);
    glp_set_col_kind(Prob, uij, GLP_BV);
    
	sprintf(buf,"u%i,%i",j,i);
	glp_set_col_name(Prob, uij+1, buf);
    glp_set_col_kind(Prob, uij+1, GLP_BV);
    
	cardinal = glp_add_rows(Prob, 2);
	
	sprintf(buf,"S%i,%i",i,j);
	glp_set_row_name(Prob, cardinal, buf);
	glp_set_row_bnds(Prob, cardinal, GLP_LO, 0, 0);
	
	constr[1] = xj;
	constr[2] = xi;
	constr[3] = uij;
	constr[4] = uij+1;
	cValues[1] = 1;
	cValues[2] = -1;
	cValues[3] = -planes[i].sep[j];
	cValues[4] = planes[i].latest - planes[j].earliest;
	
	glp_set_mat_row(Prob, cardinal, 4, constr, cValues);
	
	sprintf(buf,"S%i,%i",j,i);
	glp_set_row_name(Prob, cardinal+1, buf);
	glp_set_row_bnds(Prob, cardinal+1, GLP_LO, 0, 0);
	
	constr[1] = xi;
	constr[2] = xj;
	constr[3] = uij+1;
	constr[4] = uij;
	cValues[1] = 1;
	cValues[2] = -1;
	cValues[3] = -planes[j].sep[i];
	cValues[4] = planes[j].latest - planes[i].earliest;
	
	glp_set_mat_row(Prob, cardinal+1, 4, constr, cValues);
	
	cardinal = glp_add_rows(Prob, 1);

	sprintf(buf,"E%i,%i",i,j);
	glp_set_row_name(Prob, cardinal, buf);
	glp_set_row_bnds(Prob, cardinal, GLP_FX, 1, 0);
	
	constr[1] = uij;
	constr[2] = uij+1;
	cValues[1] = 1;
	cValues[2] = 1;
	
	glp_set_mat_row(Prob, cardinal, 2, constr, cValues);
}

/* Swap the arrival of two adjacent planes */
void swapConstraint(glp_prob * Prob, int i, int solution[], int back) {
	int t,j;
	static char firstTime = 1;
	static int u[MAXSIZE][MAXSIZE];
	char buf[AUXSIZE];
	
	if( firstTime ) {
		firstTime = 0;
		for( t = 0 ; t < n ; ++t )
		    for( j = i+1 ; j < n ; ++j ) {
            	sprintf(buf,"u%i,%i",t,j);
            	u[t][j] = glp_find_col(Prob, buf);
            	sprintf(buf,"u%i,%i",j,t);
            	u[j][t] = glp_find_col(Prob, buf);
			}
	}
	
	if( u[solution[i]][solution[i+1]] ) {
    	glp_set_col_bnds(Prob, u[solution[i]][solution[i+1]], GLP_FX, back, 0);
    	glp_set_col_bnds(Prob, u[solution[i+1]][solution[i]], GLP_FX, !back, 0);
	}
}

int compIdealT ( const void * f, const void * s ) {
	return ((struct planeOrder *)f)->ideal - ((struct planeOrder *)s)->ideal;
}
