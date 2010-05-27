#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "include\glpk.h"

#define AUX_SIZE 100
#define MAXSIZE 100
#define POSWEIGHT 10.
#define ALPHA 1.25

#define SEED 0xDEAD

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

Plane planes[MAXSIZE];

/* Creates a solution in a greedy randomized way */
int createSolution( int solution[], int time[], int pos );

/* Restriction bci = ai - bi + xi = Ti and z += ai*Ei + bi*Li */
void addBasicRestriction(glp_prob * Prob, int plane);

/* Restriction xj - xi >= Sij */
void addSeparationConstraint(glp_prob * Prob, int plane1, int plane2);

/* Restrictions xj - xi >= Sij&ij + (Li - Ej)&ji and
                xi - xj >= Sji&ji + (Lj - Ei)&ij and
				&ij + &ji = 1*/
void addOrderConstraint(glp_prob * Prob, int plane1, int plane2);

int main(void) {
	int i,j;
	
	srand(SEED);
	/* Read input */
	
	/* Order planes according to ideal time */
	
	/* Create lp instance */
	glp_prob * Prob;
	Prob = glp_create_prob();
	glp_set_prob_name(Prob, "Airplane Landing Problem");
	glp_set_obj_name(Prob, "Cost");
	
	for( i = 0 ; i < n ; ++i ) {
        addBasicRestriction(Prob,i);
	}
	
	glp_create_index(Prob);
	
	for( i = 0 ; i < n ; ++i ) {
		for( j = i+1 ; j < n ; ++j ) {
			if( planes[i].latest >= planes[j].earliest ||
			    planes[j].latest >= planes[i].earliest ) {
                addOrderConstraint(Prob,i,j);
			} else if ( planes[i].latest + planes[i].sep[j] < planes[j].earliest ) {
                addSeparationConstraint(Prob, i, j);
			} else {
                addSeparationConstraint(Prob, j, i);
			}
		}
	}
	return 0;
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
			posDiff = (pos - planes[i].pos)/POSWEIGHT;
			if( posDiff < 0 )
			    posDiff = -posDiff;
			if( planes[i].ideal - arrival > 0 )
				cost[i] = (planes[i].ideal - arrival)*planes[i].costE*(1. + posDiff);
			else
			    cost[i] = (planes[i].ideal - arrival)*planes[i].costL;
		}

	/* Find out the best option */
	double bestOption;
tryAgain:
	bestOption = DBL_MAX/(ALPHA*ALPHA);
	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && cost[i] > 0 && bestOption > cost[i] )
			bestOption = cost[i];

	/* Find out the number of options */
	int numValid;
	numValid = 0;
	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && cost[i] < bestOption * ALPHA )
			++numValid;

	/* Select an option close or equal to the best */
	int option, counter;
	option = rand()%numValid + 1;
	counter = 0;
	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && cost[i] < bestOption * ALPHA )
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

/* Restriction bci = ai - bi + xi = Ti and z += ai*Ei + bi*Li */
void addBasicRestriction(glp_prob * Prob, int plane) {
    int cardinal, constr[4], i = plane;
	double cValues[4];
	char buf[AUX_SIZE];

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
	glp_set_col_bnds(Prob, cardinal+2, GLP_LO, 0, 0);
	
    glp_add_rows(Prob, 1);
    
	sprintf(buf,"bc%i",i);
	glp_set_row_name(Prob, cardinal+3, buf);
	glp_set_row_bnds(Prob, cardinal+3, GLP_FX, planes[i].ideal, 0);
	
	constr[3] = 1 + (constr[2] = 1 + (constr[1] = cardinal));
	cValues[3] = cValues[1] = 1;
	cValues[2] = -1;
	
	glp_set_mat_row(Prob, cardinal+3, 3, constr, cValues);
}

/* Restriction xj - xi >= Sij */
void addSeparationConstraint(glp_prob * Prob, int plane1, int plane2) {
    int cardinal, constr[3], i = plane1, j = plane2;
	double cValues[3];
	char buf[AUX_SIZE];
	
    cardinal = glp_add_rows(Prob, 1);
    
	sprintf(buf,"S%i%i",i,j);
	glp_set_row_name(Prob, cardinal, buf);
	glp_set_row_bnds(Prob, cardinal, GLP_LO, planes[i].sep[j], 0);
	
	sprintf(buf,"x%i",j);
	constr[1] = glp_find_row(Prob, buf);
	
	sprintf(buf,"x%i",i);
	constr[2] = glp_find_row(Prob, buf);
	
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
	char buf[AUX_SIZE];
	int xi, xj, uij;
	
	sprintf(buf,"x%i",i);
	xi = glp_find_row(Prob, buf);
    sprintf(buf,"x%i",j);
	xj = glp_find_row(Prob, buf);
	
	uij = glp_add_cols(Prob, 2);
	
	sprintf(buf,"u%i%i",i,j);
	glp_set_col_name(Prob, uij, buf);
    glp_set_col_kind(Prob, uij, GLP_BV);
    
	sprintf(buf,"u%i%i",j,i);
	glp_set_col_name(Prob, uij+1, buf);
    glp_set_col_kind(Prob, uij+1, GLP_BV);
    
	cardinal = glp_add_rows(Prob, 2);
	
	sprintf(buf,"S%i%i",i,j);
	glp_set_row_name(Prob, cardinal, buf);
	glp_set_row_bnds(Prob, cardinal, GLP_LO, 0, 0);
	
	constr[1] = xj;
	constr[2] = xi;
	constr[3] = uij;
	constr[4] = uij+1;
	cValues[1] = 1;
	cValues[2] = -1;
	cValues[3] = -planes[i].sep[j];
	cValues[4] = planes[i].latest - planes[j].latest;
	
	glp_set_mat_row(Prob, cardinal, 4, constr, cValues);
	
	sprintf(buf,"S%i%i",j,i);
	glp_set_row_name(Prob, cardinal+1, buf);
	glp_set_row_bnds(Prob, cardinal+1, GLP_LO, 0, 0);
	
	constr[1] = xi;
	constr[2] = xj;
	constr[3] = uij+1;
	constr[4] = uij;
	cValues[1] = 1;
	cValues[2] = -1;
	cValues[3] = -planes[j].sep[i];
	cValues[4] = planes[j].latest - planes[i].latest;
	
	glp_set_mat_row(Prob, cardinal+1, 4, constr, cValues);
	
	cardinal = glp_add_rows(Prob, 1);

	sprintf(buf,"E%i%i",i,j);
	glp_set_row_name(Prob, cardinal, buf);
	glp_set_row_bnds(Prob, cardinal, GLP_FX, 1, 0);
	
	constr[1] = uij;
	constr[2] = uij+1;
	cValues[1] = 1;
	cValues[2] = -1;
	
	glp_set_mat_row(Prob, cardinal, 2, constr, cValues);
}
