#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define S 100
#define TENT 5

int main(void) {
	int i,j,n;
	char in[S], exec[S], param[S], command[S], out[S];
	int maxAlpha, randomSeed;
	FILE * file;
	
	srand(time(NULL));
	
	sprintf(exec,"Avioes.exe");
	sprintf(out,"out.txt");
	
	n = 8;
	for( i = 1 ; i <= 8 ; ++i ) {
		sprintf(in,"entradas/airland%i.txt",i);
		for( j = 0 ; j < TENT ; ++j ) {
			randomSeed = rand();
			maxAlpha = j%5 + 1;
			file = fopen(out,"a");
			fprintf(file,"%i 100 20 %i ",maxAlpha,randomSeed);
			fclose(file);
			sprintf(param,"%s %i 100 20 %i 1",out,maxAlpha,randomSeed);
			sprintf(command,"%s %s < %s",exec,param,in);
			system(command);
		}
	}
	return 0;
}
