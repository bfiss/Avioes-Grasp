#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "include\glpk.h"

#define MAXSIZE 100
#define POSWEIGHT 10.
#define ALPHA 1.25

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
	int i;
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
	    	for( j = pos - 1 ; j >= 0 ; --j )
				if( arrival < planes[solution[j]].sep[i] + time[j] )
			    	arrival = planes[solution[j]].sep[i] + time[j];
			if( arrival > planes[i].latest )
		    	return 1; /* Impossible */
			if( arrival < planes[i].earliest )
		    	arrival = planes[i].earliest;
			if( arrival >= arrivalLimit ) {
				valid[i] = 0;
				continue;
			}
			posDiff = (pos - planes[i].pos)/POSWEIGHT;
			if( posDiff < 0 )
			    posDiff = -posDiff;
			time[i] = arrival;
			if( planes[i].ideal - arrival > 0 )
				cost[i] = (planes[i].ideal - arrival)*planes[i].costE*(1. + posDiff);
			else
			    cost[i] = (planes[i].ideal - arrival)*planes[i].costL;
		}
		
	/* Find out the best option */
	double bestOption = DBL_MAX;
	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && cost[i] > 0 && bestOption > cost[i] )
			bestOption = cost[i];

	/* Find out the number of options */
	int numValid = 0;
	for( i = 0 ; i < n ; ++i )
	    if( valid[i] && cost[i] < bestOption * ALPHA )
			++numValid;

	/* Select an option close or equal to the best */
}

int main(void) {
	return 0;
}
