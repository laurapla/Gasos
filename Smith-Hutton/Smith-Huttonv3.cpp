#include <iostream>
#include <math.h>
#include<fstream>

using namespace std;

// Numerical parameters
const int N = 50;
const int M = 25;

typedef double matrix[M+2][N+2];
typedef double mface[M+1][N+1];


// FUNCTIONS
void coordinates(float dx, int N, float xvc[], float x[]);
void surface(float *yvc, int M, float Sv[]);
void volume(float dx, float dy, int N, int M, matrix& V);
void velocity(float *x, float *y, int N, int M, matrix& u, matrix& v);
void phi_inlet_outlet(float *x, float alpha, int N, double phis[]);
double max(double a, double b);
double Aperator(string method, double P);
void coefficients (int N, int M, string method, float rho0, float gamma, float dt, float Sp, float Sc, float *x, float *y, float *Sh, float *Sv, matrix V, matrix phi0, mface mflowx, mface mflowy, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap, matrix& bp);
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float* x, double* phis, double phi_boundary, float fr, float delta, int N, int M, matrix& T);
void output_matrix(int N, int M, matrix mat);
void output_file (matrix T, int N);


int main(){
	
	// DATA
	float alpha = 10; // Angle [º]
	float rho = 1; // Density
	float gamma = rho/1000000;
	float Sc = 0; // Source term = Sc+Sp*phi
	float Sp = 0;
	string method = "PLDS";
	
	float delta = 0.0000000001; // Precision of the simulation
	float fr = 1.2; // Relaxation factor
	int Time = 100;
	
	// PREVIOUS CALCULATIONS
	
	// Increments
	float dx, dy, dt;
	dx = 2.0/N;
	dy = 1.0/M;
	dt = 1;
	
	// Coordinates
	float xvc[N+1], yvc[M+1]; // Coordinates of the faces
	float x[N+2], y[M+2]; // Coordinates of the nodes
	
	xvc[0] = -1;
	x[0] = xvc[0];
	coordinates(dx, N+1, xvc, x);
	x[N+1] = 1;
	xvc[0] = -1+dx/2;
	xvc[N] = 1-dx/2;
	
	yvc[0] = 1;
	y[0] = yvc[0];
	coordinates(-dy, M+1, yvc, y);
	y[M+1] = 0;
	yvc[0] = 1-dy/2;
	yvc[M] = dy/2;
	
	// Surfaces and volumes
	float Sh[N+2], Sv[M+2];
	matrix V;
	surface(yvc, M+2, Sv);
	surface(xvc, N+2, Sh);
	volume(dx, dy, N+2, M+2, V);
	
	
	// Boundary conditions
	double phi_boundary, phis[N+1];
	phi_inlet_outlet(x, alpha, N+2, phis);
	phi_boundary = 1-tanh(alpha);
	
	
	// Mass flow on the faces
	mface mflowx, mflowy;
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+1; j++)
		{
			mflowx[j][i] = rho*Sv[j]*2*yvc[j]*(1-pow(xvc[i],2));
			mflowy[j][i] = -rho*Sh[i]*2*xvc[i]*(1-pow(yvc[j],2));
		}
	}
	
	
	// Assignation
	matrix phi, phi0;
	for(int i = 0; i<N+2; i++)
	{
		for(int j = 0; j<M+2; j++)
		{
			phi[j][i] = 1;
		}
	}
	
	
	matrix ae, aw, an, as, ap, bp;
	double resta;
	double MAX;
	float t = 0;
	
	while(t<=Time)
	{
		//New increment of time
		for(int i = 0; i<N+2; i++)
		{
			for(int j = 0; j<M+2; j++)
			{
				phi0[j][i] = phi[j][i];
			}
		}
		
		coefficients (N+2, M+2, method, rho, gamma, dt, Sp, Sc, x, y, Sh, Sv, V, phi0, mflowx, mflowy, ae, aw, an, as, ap, bp);
		Gauss_Seidel (ap, aw, ae, as, an, bp, x, phis, phi_boundary, fr, delta, N+2, M+2, phi);
		
		t = t+dt;
	}
	
	
	// SCREEN!!!!!!!!! :D
	output_matrix(N+2, M+2, phi);
		
	ofstream results;
    results.open("Resultats.dat");
    for(int k = 0; k<N+2; k++)
    {
    	results<<x[k]<<"	"<<phi[M+1][k]<<"\n";
	}
    results.close();
    
}



void coordinates(float dx, int N, float xvc[], float x[])
{
	for(int i = 0; i<N; i++)
	{
		xvc[i+1] = xvc[i]+dx;
		x[i+1] = (xvc[i+1]+xvc[i])/2;
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


void volume(float dx, float dy, int N, int M, matrix& V)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if((i==0 && j==0) || (i==0 && j==M-1) || (i==N-1 && j==0) || (i==N-1 && j==M-1))
			{
				V[j][i] = dx*dy/4;
			}
			else if(i==0 || j==0 || i==N-1 || j==M-1)
			{
				V[j][i] = dx*dy/2;
			}
			else
			{
				V[j][i] = dx*dy;
			}
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
			v[j][i] = -2*x[i]*(1-pow(y[j],2));
		}
	}
}


void phi_inlet_outlet(float *x, float alpha, int N, double phis[])
{
	for(int i = 0; i<N; i++)
	{
		if(x[i]<=0)
		{
			phis[i] = 1+tanh(alpha*(2*x[i]+1));
		}
		else
		{
			phis[i] = 0;
		}
	}
}


double max(double a, double b)
{
	if(a>b)
	{
		return a;
	}
	else
	{
		return b;
	}
}


double Aperator(string method, double P)
{
	double A;
	if(method=="CDS") // Central Differencing Scheme
	{
		A = 1-0.5*fabs(P);
	}
	else if(method=="UDS") // Upwind Differencing Scheme
	{
		A = 1;
	}
	else if(method=="HDS") // Hybrid Differencing Scheme
	{
		A = max(0,1-0.5*fabs(P));
	}
	else if(method=="PLDS") // Power Law Differencing Scheme
	{
		A = max(0,pow(1-0.1*fabs(P),5));
	}
	else if(method=="EDS") // Exponential Differencing Scheme
	{
		A = fabs(P)/(exp(fabs(P))-1);
	}
	return A;
}


void coefficients (int N, int M, string method, float rho0, float gamma, float dt, float Sp, float Sc, float *x, float *y, float *Sh, float *Sv, matrix V, matrix phi0, mface mflowx, mface mflowy, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap, matrix& bp)
{
	double De, Dw, Dn, Ds;
	double Pe, Pw, Pn, Ps;
	double Fe, Fw, Fn, Fs;
	
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(j==M-1)
			{
				ae[j][i] = 0;
				aw[j][i] = 0;
				an[j][i] = 1;
				as[j][i] = 0;
				ap[j][i] = 1;
				bp[j][i] = 0;
			}
			else
			{
				De = gamma*Sh[i]/fabs(x[i+1]-x[i]);
				Dw = gamma*Sh[i-1]/fabs(x[i]-x[i-1]);
				Dn = gamma*Sv[j-1]/fabs(y[j-1]-y[j]);
				Ds = gamma*Sv[j]/fabs(y[j]-y[j+1]);
				Fe = mflowx[j-1][i];
				Fw = mflowx[j-1][i-1];
				Fn = mflowy[j-1][i-1];
				Fs = mflowy[j][i-1];
				Pe = Fe/De;
				Pw = Fw/Dw;
				Pn = Fn/Dn;
				Ps = Fs/Ds;
				ae[j][i] = De*Aperator(method,Pe)+max(-Fe,0);
				aw[j][i] = Dw*Aperator(method,Pw)+max(Fw,0);
				an[j][i] = Dn*Aperator(method,Pn)+max(-Fn,0);
				as[j][i] = Ds*Aperator(method,Ps)+max(Fs,0);
				ap[j][i] = ae[j][i]+aw[j][i]+an[j][i]+as[j][i]+rho0*V[j][i]/dt-Sp*V[j][i];
				bp[j][i] = rho0*V[j][i]*phi0[j][i]/dt+Sc*V[j][i];
			}
		}
	}
}


// Solver (using Gauss-Seidel)
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float* x, double* phis, double phi_boundary, float fr, float delta, int N, int M, matrix& T)
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
				if(j==M-1)
				{
					if(x[i]<0)
					{
						T[j][i] = phis[i];
					}
					else
					{
						T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
				}
				else if(i==0 || j==0 || i==N-1)
				{
					T[j][i] = phi_boundary;
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


// Create an output file with the results
void output_file (matrix T, int N)
{
	ofstream results;
    results.open("Resultats.dat");
    for(int k = N/2-1; k<N; k++)
    {
    	results<<T[M-1][k]<<"\n";
	}
    results.close();
}
