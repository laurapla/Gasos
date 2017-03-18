#include <iostream>
#include <math.h>

using namespace std;

// Numerical parameters
const int N = 10;
const int M = 5;

typedef double matrix[M][N];
typedef double mface[M+1][N+1];


// FUNCTIONS
void coordinates(float dx, int N, float xvc[], float x[]);
void surface(float *yvc, int M, float Sv[]);
void volume(float *xvc, float *yvc, int N, int M, matrix& V);
void velocity(float *x, float *y, int N, int M, matrix& u, matrix& v);
void phi_inlet_outlet(float *xvc, float alpha, int N, double phis[]);
void output_matrix(int N, int M, matrix mat);


int main(){
	
	// DATA
	int alpha = 10; // Angle [º]
	float rho = 1; // Density
	float gamma = rho/10;
	float Sc = 0; // Source term = Sc+Sp*phi
	float Sp = 0;
	string method = "CDS";
	
	float delta = 1; // Precision of the simulation
	float fr = 1.2; // Relaxation factor
	
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
	
	
	// Mass flow on the faces
	mface mflowx, mflowy;
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+1; j++)
		{
			mflowx[j][i] = rho*Sv[j]*2*yvc[j]*(1-pow(xvc[i],2));
			mflowy[j][i] = rho*Sh[i]*2*xvc[i]*(1-pow(yvc[j],2));
		}
	}
	
	
	// Value on the faces
	matrix phiadim, phi;
	mface phifadim, phif;
	double phid, phic, phiu;
	float xd, xc, xu;
	float xfd, xfc, xfu;
	double xadim, xfadim;
	for(int i = 0; i<N; i++)
	{
		xc = x[i];
		xfc = xvc[i];
		if(i==0)
		{
			xu = x[i+1];
			xd = xvc[i];
			xfu = xvc[i+1];
			xfd = xvc[i-1];
		}
		else if(i==N-1)
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
				phid = phis[i];
			}
			else if(i==N-1)
			{
				phiu = phi_boundary;
				phid = phi[j][i-1];
			}
			else
			{
				phiu = phi[j][i+1];
				phid = phi[j][i-1];
			}
			phiadim[j][i] = (phic-phiu)/(phid-phiu);
			
			if(method=="CDS")
			{
				phifadim[j][i] = (xfadim-xadim)/(1-xadim)+(xfadim-1)*phiadim[j][i]/(xadim-1);
			}
			else if(method=="UDS")
			{
				phifadim[j][i] = phiadim[j][i];
			}
			else if(method=="SUDS")
			{
				phifadim[j][i] = xfadim*phiadim[j][i]/xadim;
			}
			else if(method=="QUICK")
			{
				phifadim[j][i] = xfadim+xfadim*(xfadim-1)*(phiadim[j][i]-xadim)/(xadim*(xadim-1));
			}
			
			phif[j][i] = phiu+phiadim[j][i]*(phid-phiu);
		}
	}
	
	
	// SCREEN!!!!!!!!! :D
//	output_matrix(N, M, v);
//	for(int j = 0; j<=N; j++)
//	{
//		cout<<xvc[j]<<"	";
//	}
//	for(int j = 0; j<=N; j++)
//	{
//		cout<<phis[j]<<"	";
//	}
	for(int j = 0; j<=M; j++)
	{
		for(int i = 0; i<=N; i++)
		{
			cout<<phif[j][i]<<"	";
		}
		cout<<endl;
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


void coefficients ()
{
	if(i==0 && j==0)
	{
		De = gamma*Sh[i+1]/fabs(x[i+1]-x[i]);
		Ds = gamma*Sv[j+1]/fabs(y[j]-y[j+1]);
		ae[j][i] = De-(m[j][i+1]-fabs(m[j][i+1]))/2;
		aw[j][i] = 0;
		an[j][i] = 0;
		as[j][i] = Ds+(m[j+1][i]+fabs(m[j+1][i]))/2;
		ap[j][i] = ae[j][i]+aw[j][i]+an[j][i]+as[j][i]+rho0*V[j][i]/dt-Sp*V[j][i];
		bp[j][i] = -m[j][i+1]*(phiHRS[j][i+1]-phiUDS[j][i+1])+m[j][i-1]*(phiHRS[j][i-1]-phiUDS[j][i-1])-m[j-1][i]*(phiHRS[j-1][i]-phiUDS[j-1][i])+m[j+1][i]*(phiHRS[j+1][i]-phiUDS[j+1][i])+rho0*V[j][i]*phi0[j][i]/dt+Sc*V[j][i];
	}
	else
	{
		De = gamma*Sh[i+1]/fabs(x[i+1]-x[i]);
		Dw = gamma*Sh[i]/fabs(x[i]-x[i-1]);
		Dn = gamma*Sv[j]/fabs(y[j-1]-y[j]);
		Ds = gamma*Sv[j+1]/fabs(y[j]-y[j+1]);
		ae[j][i] = De-(m[j][i+1]-fabs(m[j][i+1]))/2;
		aw[j][i] = Dw+(m[j][i-1]-fabs(m[j][i-1]))/2;
		an[j][i] = Dn-(m[j-1][i]-fabs(m[j-1][i]))/2;
		as[j][i] = Ds+(m[j+1][i]+fabs(m[j+1][i]))/2;
		ap[j][i] = ae[j][i]+aw[j][i]+an[j][i]+as[j][i]+rho0*V[j][i]/dt-Sp*V[j][i];
		bp[j][i] = -m[j][i+1]*(phiHRS[j][i+1]-phiUDS[j][i+1])+m[j][i-1]*(phiHRS[j][i-1]-phiUDS[j][i-1])-m[j-1][i]*(phiHRS[j-1][i]-phiUDS[j-1][i])+m[j+1][i]*(phiHRS[j+1][i]-phiUDS[j+1][i])+rho0*V[j][i]*phi0[j][i]/dt+Sc*V[j][i];
	}
	
}


// Solver (using Gauss-Seidel)
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float fr, float delta, int N, int M, matrix& T)
{
	double Tcalc[M][N]; // Temperature calculated in the previous iteration
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			Tcalc[j][i] = T[j][i];
		}
	}
	
	double MAX = 1; // Maximum value of the difference between T and Tcalc
	double resta = 1; // Difference between T and Tcalc
		
	while(MAX>delta)
	{			
		
		// SOLVER: Gauss-Seidel
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M; j++)
			{
				if(i==0 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==0 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==0 && j!=0 && j!=M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N-1 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N-1 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+an[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N && j!=0 && j!=M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i!=0 && i!=N-1 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i!=0 && i!=N-1 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
			}
		}
		
		// Comprovation
		MAX = 0;
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M; j++)
			{
				resta = fabs(Tcalc[j][i]-T[j][i]);
				
				if(resta>MAX)
				{
					MAX = resta;
				}
			}
		}
	
		// New assignation
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M; j++)
			{
				Tcalc[j][i] = T[j][i];
			}
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
