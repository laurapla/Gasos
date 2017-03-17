#include <iostream>
#include <math.h>

using namespace std;

// Numerical parameters
const int N = 10;
const int M = 5;

typedef double matrix[M][N];


// FUNCTIONS
void coordinates(float dx, int N, float xvc[], float x[]);
void surface(float *yvc, int M, float Sv[]);
void volume(float *xvc, float *yvc, int N, int M, matrix& V);
void velocity(float *x, float *y, int N, int M, matrix& u, matrix& v);
void phi_inlet_outlet(float *xvc, float alpha, int N, double phis[]);
void output_matrix(int N, int M, matrix mat);


int main(){
	
	// DATA
	int alpha = 10;
	float Sc = 0; // Source term = Sc+Sp*phi
	float Sp = 0;
	string method = "CDS";
	
	// PREVIOUS CALCULATIONS
	
	// Increments
	float dx, dy;
	dx = 2.0/N;
	dy = 1.0/M;
	
	// Coordinates
	float xvc[N+1], yvc[M+1]; // Coordinates of the faces
	float x[N], y[M]; // Coordinates of the nodes
	
	xvc[0] = -1;
	coordinates(dx, N, xvc, x);
	
	yvc[0] = 1;
	coordinates(-dy, M, yvc, y);
	
	// Surfaces and volumes
	float Sh[N], Sv[M];
	matrix V;
	surface(yvc, M, Sv);
	surface(xvc, N, Sh);
	volume(xvc, yvc, N, M, V);
	
	// Velocity
	matrix u, v;
	velocity(x, y, N, M, u, v);
	
	// Boundary conditions
	double phi_boundary, phis[N+1];
	phi_inlet_outlet(xvc, alpha, N, phis);
	phi_boundary = 1-tanh(alpha);
	
	
	// Value on the faces
	matrix phiadim;
	double phid, phic, phiu;
	float xd, xc, xu;
	float xfd, xfc, xfu;
	for(int i = 0; i<N; i++)
	{
		xc = x[i];
		xf = xvc[i];
		if(i==0)
		{
			xu = x[i+1];
			xd = xvc[i];
			xfu = xvc[i+1];
			xfd = xvc[i-1];
		}
		elseif(i==N-1)
		{
			xu = xvc[i+1];
			xd = x[i-1];
		}
		else
		{
			xu = x[i+1];
			xd = x[i-1];
		}
		xadim = (xc-xu)/(xd-xu);
		xfadim = (xfc-xfu)/(xfd-xfu);
		
		for(int j = 0; j<M; j++)
		{
			phic = phi[j][i];
			if(i==0)
			{
				phiu = phi[j][i+1];
				phid = phis;
			}
			elseif(i==N-1)
			{
				phiu = phi_boundary;
				phid = phi[j][i-1];
			}
			else
			{
				phiu = phi[j][i+1];
				phid = phi[j][i-1];
			}
			phiadim = (phic-phiu)/(phid-phiu);
			
			if(method=="CDS")
			{
				phifadim = 
			}
		}
	}
	
	
	// SCREEN!!!!!!!!! :D
//	output_matrix(N, M, v);
	for(int j = 0; j<=N; j++)
	{
		cout<<xvc[j]<<"	";
	}
	for(int j = 0; j<=N; j++)
	{
		cout<<phis[j]<<"	";
	}
	
}



void coordinates(float dx, int N, float xvc[], float x[])
{
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


void velocity(float *x, float *y, int N, int M, matrix& u, matrix& v)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			u[j][i] = 2*y[j]*(1-pow(x[i],2));
			v[j][i] = 2*x[i]*(1-pow(y[j],2));
		}
	}
}


void phi_inlet_outlet(float *xvc, float alpha, int N, double phis[])
{
	for(int i = 0; i<=N; i++)
	{
		if(xvc[i]<=0)
		{
			phis[i] = 1+tanh(alpha*(2*xvc[i]+1));
		}
		else
		{
			phis[i] = 0;
		}
	}
}

//
//void CDS(float *xvc, float *x, matrix phi, matrix& phif)
//{
//	double fe = fabs(xvc[j][i+1]-x[j][i])/fabs(x[j][i+1]-x[j][i]);
//	phif[j][i+1] = phi[j][i]+fe*(phi[j][i+1]-phi[j][i]);
//}
//
//
//void CDS(float *xvc, float *x, matrix m, matrix phi, matrix& phif)
//{
//	double fe;
//	if(m[j][i+1]>=0)
//	{
//		fe = 0;
//	}
//	else
//	{
//		fe = 1;
//	}
//	phif[j][i+1] = phi[j][i]+fe*(phi[j][i+1]-phi[j][i]);
//}


//void coefficients ()
//{
//	ae[j][i] = De-(m[j][i+1]-fabs(m[j][i+1]))/2;
//	aw[j][i] = Dw+(m[j][i-1]-fabs(m[j][i-1]))/2;
//	an[j][i] = Dn-(m[j-1][i]-fabs(m[j-1][i]))/2;
//	as[j][i] = Ds+(m[j+1][i]+fabs(m[j+1][i]))/2;
//	ap[j][i] = ae[j][i]+aw[j][i]+an[j][i]+as[j][i]+rho0*V[j][i]/dt-Sp*V[j][i];
//	bp[j][i] = -m[j][i+1]*(phiHRS[j][i+1]-phiUDS[j][i+1])+m[j][i-1]*(phiHRS[j][i-1]-phiUDS[j][i-1])-m[j-1][i]*(phiHRS[j-1][i]-phiUDS[j-1][i])+m[j+1][i]*(phiHRS[j+1][i]-phiUDS[j+1][i])+rho0*V[j][i]*phi0[j][i]/dt+Sc*V[j][i]
//}


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
