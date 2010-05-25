#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "include\glpk.h"

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

int main(void) {
	srand(SEED);
	return 0;
}
