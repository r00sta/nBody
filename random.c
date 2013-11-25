#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int i;
int max = 10;
int min = 0;

srand(1);

double rand_range(double min_n, double max_n)
{
    return (double)rand()/RAND_MAX * (max_n - min_n) + min_n;
}

int main()
{
	for (i=0; i<10; i++)
	{
		double r = rand_range(-1.0, 1.0);
		printf ("%f \n",r);
	}
return 0;
}
