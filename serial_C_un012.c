// seri_c_un012.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>

#define  PI 3.1415926535897932
#define T 1000000

double** calloc2D(int x, int y);
void free2D(double **u, int x, int y);

int main()
{
	int i, j, n;
	double eps0 = pow(10, -9) / (36.0 * PI);
	double mu0 = 4 * PI * pow(10, -7);
	double mm = 2, nn = 2;
	double Ly = 100.0, Lx = 100.0;
	double c = 3.0 * (pow(10, 8));
	double f = (sqrt(pow((mm * PI / Lx), 2) + pow((nn * PI / Ly), 2))) / (2 * PI * sqrt(mu0 * eps0));
	double lambda = c / f;
	double dx, dy;
	dx = lambda / 500;
	dy = dx;
	int Nx = ceil(Lx / dx);
	int Ny = ceil(Ly / dy);
	double dt = (0.9 / c) * pow(sqrt(pow(dx, -2.0) + pow(dy, -2.0)), -1.0);

	double* u_obs = (double*)calloc(T+1, sizeof(double));
	if (u_obs == NULL)
	{
		printf(" u_obs memory cannot be allocated");
		exit(-1);
	}

	double Sx = (c * dt) / dx;
	double Sy = (c * dt) / dy;
	double A = 2 * (1 - (pow(Sx, 2) + pow(Sy, 2)));
	const int xs = ceil(Nx / 4.0) + 1;	//178
	const int ys = ceil(Ny / 4.0) + 1;	//178
	const int xobs1 = ceil(Nx / 4);	//177
	const int yobs1 = ceil(Ny / 4) + 1; //178

	double** un0, ** un1, ** un2;

	un0 = calloc2D(Nx + 1, Ny + 1); // t(n-1)
	un1 = calloc2D(Nx + 1, Ny + 1); // t(n)
	un2 = calloc2D(Nx + 1, Ny + 1); // t(n+1)

	for (n = 1; n <= T - 1; n++)
	{
		for (i = 1; i <= Nx - 1; i++)
		{
			for (j = 1; j <= Ny - 1; j++)
			{
				un2[i][j] = A * un1[i][j] - un0[i][j]
					+ pow(Sx, 2) * (un1[i + 1][j] + un1[i - 1][j])
					+ pow(Sy, 2) * (un1[i][j + 1] + un1[i][j - 1]);
			}
		}
		un1[xs - 1][ys - 1] = (un1[xs - 1][ys - 1]) + cos(2.0 * PI * f * (n + 1.0) * dt);
		u_obs[n] = un1[xobs1 - 1][yobs1 - 1];


		int ii, jj;
		for (ii = 1; ii <= Nx-1; ii++)
		{
			for (jj = 1; jj <= Ny-1; jj++)
			{
				un0[ii][jj] = un1[ii][jj];
				un1[ii][jj] = un2[ii][jj];
			}
		}
	}

	for (int n = 0; n < T; n++)
	{
		printf("uobs[%d] = %.16lf\n ", n, u_obs[n]);
	}

    free2D(uu, Nx, Ny);

	return 0;
}


double** calloc2D(int x, int y)
{

	double** u = (double**)calloc(x, sizeof(double));

	if (u == NULL)
	{
		printf(" u memory cannot be allocated inside function memalloc2D");
		exit(-1);
	}

	int i;

	for (i = 0; i < x; ++i)
	{
		u[i] = (double*)calloc(y, sizeof(double));
		//	    printf("i=%d\n",i);
		if (u[i] == NULL)
		{
			printf(" u memory cannot be allocated");
			printf("i=%d\n", i);
			exit(-1);
		}
	}
	return u;
}

void free2D(double** u, int x, int y)
{

	int i;

	for (i = 0; i < x; ++i)
	{
		free(u[i]);
	}
	free(u);
}






