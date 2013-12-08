#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
	int i;
	int j;
	int n = 10;
	printf("\nNumber of particles: %i \n\n",n);
	
	for (i=0; i<(n-1); i++)
	{
		for (j=(i+1); j<(n); j++)
		{
			printf("%i,%i\n",i,j);
		}
	}
	return 0;
}
