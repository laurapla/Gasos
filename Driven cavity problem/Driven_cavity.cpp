#include<iostream>
#include<math.h>

using namespace std;

// Numerical parameters
const int N = 10;
const int M = 10;

typedef double matrix[M][N];

void coordinates(float dx, int N, float xvc[], float x[]);
void surface(float *yvc, int M, float Sv[]);
void volume(float *xvc, float *yvc, int N, int M, matrix& V);
void output_matrix(int N, int M, matrix mat);

int main()
{
	float L = 1; // Length of the cavity
	float rho = 1; // Density
	float mu = 0.00001; // Viscosity
	float uref = 0.001; // Reference velocity
	
	float delta = 0.001; // Precision of the simulation
	float fr = 1.1; // Relaxation factor
	
	// Coordinates
	float xvc[N+1], yvc[M+1], x[N], y[M];
	coordinates(L, N, xvc, x);
	coordinates(L, M, yvc, y);
	
	// Surface and volume
	float Sh[N], Sv[M];
	surface(xvc, N, Sh); // Horizontal surface
	surface(yvc, M, Sv); // Vertical surface
	matrix V;
	volume(xvc, yvc, N, M, V); // Volume
	
	
	// Properties that are going to be calculated
	matrix p, u, v;
	
	
	resta = 1;
	
	while(resta>delta)
	{
		
	}
	
	
//	for(int i = 0; i<M+1; i++)
//	{
//		cout<<yvc[i]<<endl;
//	}
	output_matrix(N, M, V);
	
	return 0;
}



void coordinates(float L, int N, float xvc[], float x[])
{
	float dx = L/N;
	xvc[0] = 0;
	for(int i = 0; i<N; i++)
	{
		xvc[i+1] = xvc[i]+dx;
		x[i] = (xvc[i+1]+xvc[i])/2;
	}
}


void surface(float *yvc, int M, float Sv[])
{
	for(int j = 0; j<M; j++)
	{
		Sv[j] = fabs(yvc[j]-yvc[j+1]);
		if(j==M-1)
		{
			Sv[j] = Sv[j-1];
		}
	}
}


void volume(float *xvc, float *yvc, int N, int M, matrix& V)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			V[j][i] = fabs(xvc[i+1]-xvc[i])*fabs(yvc[j]-yvc[j+1]);
		}
	}
}


void output_matrix(int N, int M, matrix mat)
{
	for(int j = 0; j<M; j++)
	{
		for(int i = 0; i<N; i++)
		{
			cout<<mat[j][i]<<"	";
		}
		cout<<endl;
	}
}
