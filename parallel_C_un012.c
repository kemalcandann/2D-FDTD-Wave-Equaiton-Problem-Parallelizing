#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <mpi.h>

#define MASTER 0
#define WORKER 1
#define  PI 3.1415926535897932
#define T 1000000


int main(int argc, char* argv[])
{
	int numofprocs;
	int id;
	double t1, t2;
	double eps0 = pow(10.0, -9.0) / (36.0 * PI);
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
	double Sx = (c * dt) / dx;
	double Sy = (c * dt) / dy;
	double A = 2 * (1 - (pow(Sx, 2) + pow(Sy, 2)));
	const int xs = ceil(Nx / 4.0) + 1;	  //178
	const int ys = ceil(Ny / 4.0) + 1;	  //178	
	const int xobs1 = ceil(Nx / 4);	      //177
	const int yobs1 = ceil(Ny / 4) + 1;   //178

	double* u_obs = (double*)calloc(T + 1, sizeof(double));

	MPI_Status status;
	MPI_Init(&argc, &argv);
	t1 = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &numofprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	MPI_Comm comm;
	int dim[2], period[2], reorder;
	int coord[2];

	//Topology
	int h = 4; 
	int g = 2;

	if (h != 1 && g == 1 && numofprocs != 1)
	{
		dim[0] = h;
		period[0] = 1; period[1] = 0;
		reorder = 1;
		MPI_Cart_create(MPI_COMM_WORLD, 1, dim, period, reorder, &comm);

		//(Nx - 1) -> resimden çerçeve çıkartılıp resme odaklanılır.(Nx+1 olmalıydı ama bi alttan bi üstteen çerçeve çıkartıldı.) (Nx + 1  == 708)
		int rx = (Nx - 1) % h;
		int nx = 0;
		if (!rx)
			nx = ((Nx - 1) / h) + 2; //+2 resime sınır koşulu + ghost eklenir.
		else
			nx = ((Nx - 1 - rx) / h) + 2;

		int nx_0 = nx + rx;

		int ny = Ny + 1;

		double* un0;
		double* un1;
		double* un2;
		int oned_a;
		int oned_b;

		if (id == 0)
		{
			nx = nx_0;
		}

		un0 = (double*)calloc(nx * ny, sizeof(double));
		un1 = (double*)calloc(nx * ny, sizeof(double));
		un2 = (double*)calloc(nx * ny, sizeof(double));
		oned_a = (nx - 2) * ny;
		oned_b = (nx - 1) * ny;


		for (int n = 1; n <= T - 1; n++) // 
		{
			if (id == 0)
			{
				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, 1, 1, comm); // send data to 1
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, 1, 2, comm); // send data to 1

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, 1, 3, comm, &status); // receive data from 1
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, 1, 4, comm, &status); // receive data from 1
			}

			else if (id == numofprocs - 1)
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, numofprocs - 2, 3, comm);
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, numofprocs - 2, 4, comm);

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, numofprocs - 2, 1, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, numofprocs - 2, 2, comm, &status);
			}

			else
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, id - 1, 3, comm); // send data to 0
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, id - 1, 4, comm); // send data to 0

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, id - 1, 1, comm, &status); // receive data from 0	
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, id - 1, 2, comm, &status); // receive data from 0	

				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, id + 1, 1, comm);
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, id + 1, 2, comm);

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, id + 1, 3, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, id + 1, 4, comm, &status);
			}

				for (int i = 1; i <= nx - 2; i++)
				{
					for (int j = 1; j <= ny - 2; j++)
					{
						int oned_c = (i * ny) + j;
						un2[oned_c] = A * un1[oned_c] - un0[oned_c]
							+ pow(Sx, 2) * (un1[oned_c + ny] + un1[oned_c - ny])
							+ pow(Sy, 2) * (un1[oned_c + 1] + un1[oned_c - 1]);
					}
				}
			

				if (h == 4)
				{
					if (id == 0)
					{
						int xsloc = (xs - 1);
						int ysloc = (ys - 1);
						int oned_d = (xsloc * ny) + ysloc;
						un1[oned_d] = (un1[oned_d]) + cos(2.0 * PI * f * (n + 1.0) * dt);

						int xobs1loc = (xobs1 - 1);
						int yobs1loc = (yobs1 - 1);
						int oned_e = ((xobs1loc * ny) + yobs1loc);
						u_obs[n] = un1[oned_e];
					}
				}
				else if (h == 8)
				{
					if (id == 1)
					{
						int xsloc = (xs + 1 - nx_0);
						int ysloc = (ys - 1);
						int oned_d = (xsloc * ny) + ysloc;
						un1[oned_d] = (un1[oned_d]) + cos(2.0 * PI * f * (n + 1.0) * dt);

						int xobs1loc = (xobs1 + 1 - nx_0);
						int yobs1loc = (yobs1 - 1);
						int oned_e = ((xobs1loc * ny) + yobs1loc);
						u_obs[n] = un1[oned_e];
					}
				}
				

			// UPDATE the time sequences.
			{
				for (int ii = 1; ii <= nx - 2; ii++)
				{
					for (int jj = 1; jj <= ny - 2; jj++)
					{
						int k = (ii * ny) + jj;
						un0[k] = un1[k];
						un1[k] = un2[k];
					}
				}
			}

		} // end of T loop
		if (h == 4)
		{
			if (id == 0)
			{
				for (int n = T - 1; n < T; n++)
				{
					printf("u_obs[ %d ] = %.16lf\n", n, u_obs[n]);
				}
			}
		}
		else if (h == 8)
		{
			if (id == 1)
			{
				for (int n = T - 1; n < T; n++)
				{
					printf("u_obs[ %d ] = %.16lf\n", n, u_obs[n]);
				}
			}
		}

		t2 = MPI_Wtime();
		printf(" - %d - total time is %lf", id, t2 - t1);
		MPI_Finalize();
		return 0;
	}

	if (h == 1 && g != 1 && numofprocs != 1)
	{
		dim[0] = g;
		period[0] = 1; period[1] = 0;
		reorder = 1;
		MPI_Cart_create(MPI_COMM_WORLD, 1, dim, period, reorder, &comm);

		int nx = Nx + 1;

		int ry = (Ny - 1) % g;
		int ny = 0;
		if (!ry)
			ny = ((Ny - 1) / g) + 2;
		else
			ny = ((Ny - 1 - ry) / g) + 2;

		int ny_0 = ny + ry;

		double* un0;
		double* un1;
		double* un2;
		int oned_a;
		int oned_b;


		if (id == 0)
		{
			ny = ny_0;
		}

		un0 = (double*)calloc(nx * ny, sizeof(double));
		un1 = (double*)calloc(nx * ny, sizeof(double));
		un2 = (double*)calloc(nx * ny, sizeof(double));

		MPI_Datatype columntype;
		MPI_Type_vector(nx, 1, ny, MPI_DOUBLE, &columntype);
		MPI_Type_commit(&columntype);

		for (int n = 1; n <= T - 1; n++)
		{
			if (id == 0)
			{
				MPI_Send(&un0[ny - 2], 1, columntype, 1, 1, comm);   //[a][b] == [a*ny+b]
				MPI_Send(&un1[ny - 2], 1, columntype, 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, 1, 4, comm, &status);
			}

			else if (id == numofprocs - 1)
			{
				MPI_Send(&un0[1], 1, columntype, numofprocs - 2, 3, comm);
				MPI_Send(&un1[1], 1, columntype, numofprocs - 2, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, numofprocs - 2, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, numofprocs - 2, 2, comm, &status);
			}

			else
			{
				MPI_Send(&un0[1], 1, columntype, id - 1, 3, comm); // send data to 0
				MPI_Send(&un1[1], 1, columntype, id - 1, 4, comm); // send data to 0

				MPI_Recv(&un0[0], 1, columntype, id - 1, 1, comm, &status); // receive data from 0	
				MPI_Recv(&un1[0], 1, columntype, id - 1, 2, comm, &status); // receive data from 0	

				MPI_Send(&un0[ny - 2], 1, columntype, id + 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, id + 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, id + 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, id + 1, 4, comm, &status);
			}


			for (int i = 1; i <= nx - 2; i++)
			{
				for (int j = 1; j <= ny - 2; j++)
				{
					int oned_c = (i * ny) + j;
					un2[oned_c] = A * un1[oned_c] - un0[oned_c]
						+ pow(Sx, 2) * (un1[oned_c + ny] + un1[oned_c - ny])
						+ pow(Sy, 2) * (un1[oned_c + 1] + un1[oned_c - 1]);
				}
			}
			

			if (g == 4)
			{
				if (id == 0)
				{
					int xsloc = (xs - 1);
					int ysloc = (ys - 1);
					int oned_d = (xsloc * ny) + ysloc;
					un1[oned_d] = (un1[oned_d]) + cos(2.0 * PI * f * (n + 1.0) * dt);

					int xobs1loc = (xobs1 - 1);
					int yobs1loc = (yobs1 - 1);
					int oned_e = ((xobs1loc * ny) + yobs1loc);
					u_obs[n] = un1[oned_e];
				}
			}
			else if (g == 8)
			{
				if (id == 1)
				{
					int xsloc = (xs - 1);
					int ysloc = (ys + 1 - ny_0);
					int oned_d = (xsloc * ny) + ysloc;
					un1[oned_d] = (un1[oned_d]) + cos(2.0 * PI * f * (n + 1.0) * dt);

					int xobs1loc = (xobs1 - 1);
					int yobs1loc = (yobs1 + 1 - ny_0);
					int oned_e = ((xobs1loc * ny) + yobs1loc);
					u_obs[n] = un1[oned_e];
				}
			}

			// UPDATE the time sequences.
				for (int ii = 1; ii <= nx - 2; ii++)
				{
					for (int jj = 1; jj <= ny - 2; jj++)
					{
						int k = (ii * ny) + jj;
						un0[k] = un1[k];
						un1[k] = un2[k];
					}
				}
			
		}

		if (g == 4)
		{
			if (id == 0)
			{
				for (int n = T - 1; n < T; n++)
				{
					printf("u_obs[ %d ] = %.16lf\n", n, u_obs[n]);
				}
			}
		}
		else if (g == 8)
		{
			if (id == 1)
			{
				for (int n = T - 1; n < T; n++)
				{
					printf("u_obs[ %d ] = %.16lf\n", n, u_obs[n]);
				}
			}
		}

		t2 = MPI_Wtime();
		printf(" - %d - total time is %lf", id, t2 - t1);
		MPI_Finalize();
		return 0;
	}

	if (h == 2 && g == 2 && numofprocs != 1)
	{
		dim[0] = h; dim[1] = g;
		period[0] = 1; period[1] = 0;
		reorder = 1;
		MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

		//(Nx - 1) -> resimden çerçeve çıkartılıp resme odaklanılır.(Nx+1 olmalıydı ama bi alttan bi üstteen çerçeve çıkartıldı.)
		int rx = (Nx - 1) % h;
		int nx = 0;
		if (!rx)
			nx = ((Nx - 1) / h) + 2; //+2 resime sınır koşulu + ghost eklenir.
		else
			nx = ((Nx - 1 - rx) / h) + 2;

		int nx_0 = nx + rx;


		int ry = (Ny - 1) % g;
		int ny = 0;
		if (!ry)
			ny = ((Ny - 1) / g) + 2;
		else
			ny = ((Ny - 1 - ry) / g) + 2;

		int ny_0 = ny + ry;

		double* un0;
		double* un1;
		double* un2;
		int oned_a;
		int oned_b;

		if (id == 0)
		{
			nx = nx_0;
			ny = ny_0;
		}
		else if (id == 1)
		{
			nx = nx_0;
		}
		else if (id == 2)
		{
			ny = ny_0;
		}
		else
		{

		}

		un0 = (double*)calloc(nx * ny, sizeof(double));
		un1 = (double*)calloc(nx * ny, sizeof(double));
		un2 = (double*)calloc(nx * ny, sizeof(double));
		oned_a = (nx - 2) * ny;
		oned_b = (nx - 1) * ny;

		MPI_Datatype columntype;
		MPI_Type_vector(nx, 1, ny, MPI_DOUBLE, &columntype);
		MPI_Type_commit(&columntype);

		for (int n = 1; n <= T - 1; n++)
		{
			if (id == 0)
			{
				MPI_Send(&un0[ny - 2], 1, columntype, 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, 1, 4, comm, &status);


				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, 2, 5, comm); //un0[nx-2][0] == un0[(nx -2)*ny + 0]
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, 2, 6, comm);

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, 2, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, 2, 8, comm, &status);
			}

			else if (id == 1)
			{
				MPI_Send(&un0[1], 1, columntype, 0, 3, comm);
				MPI_Send(&un1[1], 1, columntype, 0, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, 0, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, 0, 2, comm, &status);


				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, 3, 9, comm); // send data to 1
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, 3, 10, comm); // send data to 1

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, 3, 11, comm, &status); // receive data from 1
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, 3, 12, comm, &status); // receive data from 1
			}

			else if (id == 2)
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, 0, 7, comm); // send data to 0  //un0[1][0] == un0[(1*ny) + 0]
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, 0, 8, comm); // send data to 0

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, 0, 5, comm, &status); // receive data from 0	
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, 0, 6, comm, &status); // receive data from 0


				MPI_Send(&un0[ny - 2], 1, columntype, 3, 13, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, 3, 14, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, 3, 15, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, 3, 16, comm, &status);

			}
			else
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, 1, 11, comm); // send data to 0
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, 1, 12, comm); // send data to 0

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, 1, 9, comm, &status); // receive data from 0	
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, 1, 10, comm, &status); // receive data from 0


				MPI_Send(&un0[1], 1, columntype, 2, 15, comm);
				MPI_Send(&un1[1], 1, columntype, 2, 16, comm);

				MPI_Recv(&un0[0], 1, columntype, 2, 13, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, 2, 14, comm, &status);
			}


			for (int i = 1; i <= nx - 2; i++)
			{
				for (int j = 1; j <= ny - 2; j++)
				{
					int oned_c = (i * ny) + j;
					un2[oned_c] = A * un1[oned_c] - un0[oned_c]
						+ pow(Sx, 2) * (un1[oned_c + ny] + un1[oned_c - ny])
						+ pow(Sy, 2) * (un1[oned_c + 1] + un1[oned_c - 1]);
					//printf("un0[10] = %lf\n", un0[10]);
					//printf("un1[10] = %lf\n", un1[10]);
					//printf("un2[10] = %lf\n", un2[10]);
				}
			}



			if (id == 0)
			{
				int xsloc = (xs - 1);
				int ysloc = (ys - 1);
				int oned_d = (xsloc * ny) + ysloc;
				un1[oned_d] = (un1[oned_d]) + cos(2.0 * PI * f * (n + 1.0) * dt);

				int xobs1loc = (xobs1 - 1);
				int yobs1loc = (yobs1 - 1);
				int oned_e = ((xobs1loc * ny) + yobs1loc);
				u_obs[n] = un1[oned_e];
			}

			// UPDATE the time sequences.
			for (int ii = 1; ii <= nx - 2; ii++)
			{
				for (int jj = 1; jj <= ny - 2; jj++)
				{
					int k = (ii * ny) + jj;
					un0[k] = un1[k];
					un1[k] = un2[k];
				}
			}
		}

		if (id == 0)
		{
			for (int n = T - 1; n < T; n++)
			{
				printf("u_obs[ %d ] = %.16lf\n", n, u_obs[n]);
			}
		}

		t2 = MPI_Wtime();
		printf(" - %d - total time is %lf", id, t2 - t1);
		MPI_Finalize();
		return 0;
	}

	if (h == 2 && g != 2 && numofprocs != 1)
	{
		dim[0] = h; dim[1] = g;
		period[0] = 1; period[1] = 0;
		reorder = 1;
		MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

		//(Nx - 1) -> resimden çerçeve çıkartılıp resme odaklanılır.(Nx+1 olmalıydı ama bi alttan bi üstteen çerçeve çıkartıldı.)
		int rx = (Nx - 1) % h;
		int nx = 0;
		if (!rx)
			nx = ((Nx - 1) / h) + 2; //+2 resime sınır koşulu + ghost eklenir.
		else
			nx = ((Nx - 1 - rx) / h) + 2;

		int nx_0 = nx + rx;


		int ry = (Ny - 1) % g;
		int ny = 0;
		if (!ry)
			ny = ((Ny - 1) / g) + 2;
		else
			ny = ((Ny - 1 - ry) / g) + 2;

		int ny_0 = ny + ry;

		double* un0;
		double* un1;
		double* un2;
		int oned_a;
		int oned_b;

		if ((id == 1) || (id == 2) || (id == 3))
		{
			nx = nx_0;
		}
		else if ((id == 5) || (id == 6) || (id == 7))
		{

		}
		else if (id == 0)
		{
			nx = nx_0;
			ny = ny_0;
		}
		else
		{
			ny = ny_0;
		}

		un0 = (double*)calloc(nx * ny, sizeof(double));
		un1 = (double*)calloc(nx * ny, sizeof(double));
		un2 = (double*)calloc(nx * ny, sizeof(double));
		oned_a = (nx - 2) * ny;
		oned_b = (nx - 1) * ny;

		MPI_Datatype columntype;
		MPI_Type_vector(nx, 1, ny, MPI_DOUBLE, &columntype);
		MPI_Type_commit(&columntype);

		for (int n = 1; n <= T - 1; n++)
		{
			if (id == 0)
			{
				MPI_Send(&un0[ny - 2], 1, columntype, 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, 1, 4, comm, &status);


				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, g, 5, comm); //un0[nx-2][0] == un0[(nx -2)*ny + 0]
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, g, 6, comm);

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, g, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, g, 8, comm, &status);
			}

			else if (id == g - 1)
			{
				MPI_Send(&un0[1], 1, columntype, g - 2, 3, comm);
				MPI_Send(&un1[1], 1, columntype, g - 2, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, g - 2, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, g - 2, 2, comm, &status);


				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, numofprocs - 1, 5, comm);
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, numofprocs - 1, 6, comm);

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, numofprocs - 1, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, numofprocs - 1, 8, comm, &status);
			}

			else if (id == g)
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, 0, 7, comm);  //un0[1][0] == un0[(1*ny) + 0]
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, 0, 8, comm);

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, 0, 5, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, 0, 6, comm, &status);


				MPI_Send(&un0[ny - 2], 1, columntype, g + 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, g + 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, g + 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, g + 1, 4, comm, &status);
			}
			else if (id == numofprocs - 1)
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, g - 1, 7, comm);
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, g - 1, 8, comm);

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, g - 1, 5, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, g - 1, 6, comm, &status);


				MPI_Send(&un0[1], 1, columntype, numofprocs - 2, 3, comm);
				MPI_Send(&un1[1], 1, columntype, numofprocs - 2, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, numofprocs - 2, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, numofprocs - 2, 2, comm, &status);
			}

			else if (id > 0 && id < g - 1)
			{
				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, g + id, 5, comm);
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, g + id, 6, comm);

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, g + id, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, g + id, 8, comm, &status);


				MPI_Send(&un0[ny - 2], 1, columntype, id + 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, id + 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, id + 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, id + 1, 4, comm, &status);

				MPI_Send(&un0[1], 1, columntype, id - 1, 3, comm);
				MPI_Send(&un1[1], 1, columntype, id - 1, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, id - 1, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, id - 1, 2, comm, &status);
			}

			else if (id > g && id < numofprocs - 1)
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, id - g, 7, comm);
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, id - g, 8, comm);

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, id - g, 5, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, id - g, 6, comm, &status);


				MPI_Send(&un0[ny - 2], 1, columntype, id + 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, id + 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, id + 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, id + 1, 4, comm, &status);

				MPI_Send(&un0[1], 1, columntype, id - 1, 3, comm);
				MPI_Send(&un1[1], 1, columntype, id - 1, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, id - 1, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, id - 1, 2, comm, &status);
			}


			for (int i = 1; i <= nx - 2; i++)
			{
				for (int j = 1; j <= ny - 2; j++)
				{
					int oned_c = (i * ny) + j;
					un2[oned_c] = A * un1[oned_c] - un0[oned_c]
						+ pow(Sx, 2) * (un1[oned_c + ny] + un1[oned_c - ny])
						+ pow(Sy, 2) * (un1[oned_c + 1] + un1[oned_c - 1]);
				}
			}



			if (id == 0)
			{
				int xsloc = (xs - 1);
				int ysloc = (ys - 1);
				int oned_d = (xsloc * ny) + ysloc;
				un1[oned_d] = (un1[oned_d]) + cos(2.0 * PI * f * (n + 1.0) * dt);

				int xobs1loc = (xobs1 - 1);
				int yobs1loc = (yobs1 - 1);
				int oned_e = ((xobs1loc * ny) + yobs1loc);
				u_obs[n] = un1[oned_e];
			}
			
			// UPDATE the time sequences.
			for (int ii = 1; ii <= nx - 2; ii++)
			{
				for (int jj = 1; jj <= ny - 2; jj++)
				{
					int k = (ii * ny) + jj;
					un0[k] = un1[k];
					un1[k] = un2[k];
				}
			}
		}

		if (id == 0)
		{
			for (int n = T - 1; n < T; n++)
			{
				printf("u_obs[ %d ] = %.16lf\n", n, u_obs[n]);
			}
		}

		t2 = MPI_Wtime();
		printf(" - %d - total time is %lf", id, t2 - t1);
		MPI_Finalize();
		return 0;
	}

	if (h != 2 && g == 2 && numofprocs != 1)
	{
		dim[0] = h; dim[1] = g;
		period[0] = 1; period[1] = 0;
		reorder = 1;
		MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

		//(Nx - 1) -> resimden çerçeve çıkartılıp resme odaklanılır.(Nx+1 olmalıydı ama bi alttan bi üstteen çerçeve çıkartıldı.)
		int rx = (Nx - 1) % h;
		int nx = 0;
		if (!rx)
			nx = ((Nx - 1) / h) + 2; //+2 resime sınır koşulu + ghost eklenir.
		else
			nx = ((Nx - 1 - rx) / h) + 2; 

		int nx_0 = nx + rx;


		int ry = (Ny - 1) % g;
		int ny = 0;
		if (!ry)
			ny = ((Ny - 1) / g) + 2;
		else
			ny = ((Ny - 1 - ry) / g) + 2;

		int ny_0 = ny + ry;

		double* un0;
		double* un1;
		double* un2;
		int oned_a;
		int oned_b;

		if ((id == 2) || (id == 4) || (id == 6))
		{
			ny = ny_0;
		}
		else if ((id == 3) || (id == 5) || (id == 7))
		{
			
		}
		else if (id == 0)
		{
			nx = nx_0;
			ny = ny_0;
		}
		else
		{
			nx = nx_0;
		}

			un0 = (double*)calloc(nx * ny, sizeof(double));
			un1 = (double*)calloc(nx * ny, sizeof(double));
			un2 = (double*)calloc(nx * ny, sizeof(double));
			oned_a = (nx - 2) * ny;
			oned_b = (nx - 1) * ny;

		MPI_Datatype columntype;
		MPI_Type_vector(nx, 1, ny, MPI_DOUBLE, &columntype);
		MPI_Type_commit(&columntype);

		for (int n = 1; n <= T - 1; n++)
		{
			if (id == 0)
			{
				MPI_Send(&un0[ny - 2], 1, columntype, 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, 1, 2, comm);

				MPI_Recv(&un0[ny - 1], 1, columntype, 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, 1, 4, comm, &status);


				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, g, 5, comm); //un0[nx-2][0] == un0[(nx -2)*ny + 0]
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, g, 6, comm);
										 
				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, g, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, g, 8, comm, &status);
			}

			else if (id == 1)
			{
				MPI_Send(&un0[1], 1, columntype, 0, 3, comm);
				MPI_Send(&un1[1], 1, columntype, 0, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, 0, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, 0, 2, comm, &status);


				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, id + g, 5, comm);
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, id + g, 6, comm);

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, id + g, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, id + g, 8, comm, &status);
			}

			else if (id == numofprocs - 2)
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, id - g, 7, comm);  //un0[1][0] == un0[(1*ny) + 0]
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, id - g, 8, comm);

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, id - g, 5, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, id - g, 6, comm, &status);


				MPI_Send(&un0[ny - 2], 1, columntype, id + 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, id + 1, 2, comm);
								
				MPI_Recv(&un0[ny - 1], 1, columntype, id + 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, id + 1, 4, comm, &status);

			}
			else if (id == numofprocs - 1)
			{
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, id - g, 7, comm);
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, id - g, 8, comm);

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, id - g, 5, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, id - g, 6, comm, &status);


				MPI_Send(&un0[1], 1, columntype, id - 1, 3, comm);
				MPI_Send(&un1[1], 1, columntype, id - 1, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, id - 1, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, id - 1, 2, comm, &status);
			}

			else if (id > 1 && id < numofprocs - 2 && id % 2 == 0)
			{
				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, g + id, 5, comm);
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, g + id, 6, comm);
										 
				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, g + id, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, g + id, 8, comm, &status);
										 
				MPI_Send(&un0[ny], ny, MPI_DOUBLE, id - g, 7, comm);
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, id - g, 8, comm);
										 
				MPI_Recv(&un0[0], ny, MPI_DOUBLE, id - g, 5, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, id - g, 6, comm, &status);
										 

				MPI_Send(&un0[ny - 2], 1, columntype, id + 1, 1, comm);
				MPI_Send(&un1[ny - 2], 1, columntype, id + 1, 2, comm);
								
				MPI_Recv(&un0[ny - 1], 1, columntype, id + 1, 3, comm, &status);
				MPI_Recv(&un1[ny - 1], 1, columntype, id + 1, 4, comm, &status);

			}

			else if (id > 1 && id < numofprocs - 2 && id % 2 != 0)
			{
				MPI_Send(&un0[oned_a], ny, MPI_DOUBLE, g + id, 5, comm);
				MPI_Send(&un1[oned_a], ny, MPI_DOUBLE, g + id, 6, comm);

				MPI_Recv(&un0[oned_b], ny, MPI_DOUBLE, g + id, 7, comm, &status);
				MPI_Recv(&un1[oned_b], ny, MPI_DOUBLE, g + id, 8, comm, &status);

				MPI_Send(&un0[ny], ny, MPI_DOUBLE, id - g, 7, comm);
				MPI_Send(&un1[ny], ny, MPI_DOUBLE, id - g, 8, comm);

				MPI_Recv(&un0[0], ny, MPI_DOUBLE, id - g, 5, comm, &status);
				MPI_Recv(&un1[0], ny, MPI_DOUBLE, id - g, 6, comm, &status);


				MPI_Send(&un0[1], 1, columntype, id - 1, 3, comm);
				MPI_Send(&un1[1], 1, columntype, id - 1, 4, comm);

				MPI_Recv(&un0[0], 1, columntype, id - 1, 1, comm, &status);
				MPI_Recv(&un1[0], 1, columntype, id - 1, 2, comm, &status);
			}

				for (int i = 1; i <= nx - 2; i++)
				{
					for (int j = 1; j <= ny - 2; j++)
					{
						int oned_c = (i * ny) + j;
						un2[oned_c] = A * un1[oned_c] - un0[oned_c]  
							+ pow(Sx, 2) * (un1[oned_c + ny] + un1[oned_c - ny]) // un1[i][j] == un1[oned_c] 
							+ pow(Sy, 2) * (un1[oned_c + 1] + un1[oned_c - 1]);
					}
				}
			

			if (id == 0)
			{
				int xsloc = (xs - 1);
				int ysloc = (ys - 1);
				int oned_d = (xsloc * ny) + ysloc;
				un1[oned_d] = (un1[oned_d]) + cos(2.0 * PI * f * (n + 1.0) * dt);
	
				int xobs1loc = (xobs1 - 1);
				int yobs1loc = (yobs1 - 1);
				int oned_e = ((xobs1loc * ny) + yobs1loc);
				u_obs[n] = un1[oned_e];
			} 

			// UPDATE the time sequences.

				for (int ii = 1; ii <= nx - 2; ii++)
				{
					for (int jj = 1; jj <= ny - 2; jj++)
					{
						int k = (ii * ny) + jj;
						un0[k] = un1[k];
						un1[k] = un2[k];
					}
				}
			}

		if (id == 0)
		{
			for (int n = T - 1; n < T; n++)
			{
				printf("u_obs[ %d ] = %.16lf\n", n, u_obs[n]);
			}
		}

	
		t2 = MPI_Wtime();
		printf(" - %d - total time is %lf", id, t2 - t1);
		MPI_Finalize();
		return 0;
	}

}
