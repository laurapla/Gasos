#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

// Numerical parameters
const int N = 100;
const int M = 100;

typedef double matrix[M+2][N+2];
typedef double staggx[M+2][N+1];
typedef double staggy[M+1][N+2];

void coordinates(float L, int N, float xvc[], float x[]);
void surface(float *yvc, int M, float Sv[]);
void volume(float *xvc, float *yvc, int N, int M, matrix& V);
void constant_coefficients(int N, int M, float *x, float *y, float *Sv, float *Sh, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap);
double convective_term (string method, float xf, float x2, float x3, double u2, double u3);
void intermediate_velocities (int N, int M, float rho, float mu, string method, float delta, float dt, float* x, float* y, float *xvc, float* yvc, float* Sh, float* Sv, matrix V, staggx u0, staggy v0, staggx Ru0, staggy Rv0, staggx &Ru, staggy &Rv, staggx &up, staggy &vp);
void bp_coefficient (int N, int M, float rho, float dt, float* Sh, float* Sv, staggx up, staggy vp, matrix bp);
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float fr, float delta, int N, int M, matrix& T);
void velocities (int N, int M, float rho, float dt, float uref, float* x, float* y, matrix p, staggx up, staggy vp, staggx &u, staggy &v);
double min(double a, double b);
double max(double a, double b);
double time_step (float dtd, float* x, float* y, staggx u, staggy v);
void output_files (float* x, float* y, staggx u, staggy v);


int main()
{
	int Re = 100; // Reynolds number
	float L = 1; // Length of the cavity
	float rho = 1; // Density
	float uref = 1; // Reference velocity
	float mu = rho*uref*L/Re; // Viscosity
	
	
	string method = "CDS";
	float delta = 1e-2; // Precision of the simulation (error/dt)
	float fr = 1; // Relaxation factor
	
	// Coordinates
	float xvc[N+1], yvc[M+1], x[N+2], y[M+2];
	coordinates(L, N, xvc, x);
	coordinates(L, M, yvc, y);
	
	// Surfaces
	float Sh[N+2], Sv[M+2];
	matrix V;
	surface(xvc, N+2, Sh); // Horizontal surface
	surface(yvc, M+2, Sv); // Vertical surface
	volume(xvc, yvc, N+2, M+2, V); // Volume
	
	
	// Properties that are going to be calculated
	matrix p; // Values in the nodes (pressure)
	staggx u, u0, Ru, Ru0; // Values in the points given by the staggered meshes (velocities)
	staggy v, v0, Rv, Rv0;
	
	// Inicialization
	for(int j = 0; j<M+2; j++)
	{
		for(int i = 0; i<N+1; i++)
		{
			if(j==M+1 && i!=0 && i!=N)
			{
				u0[j][i] = uref; // Horizontal velocity at n
			}
			else
			{
				u0[j][i] = 0; // Horizontal velocity at n
			}
			Ru0[j][i] = 0; // R (horizontal) at n-1
		}
	}
	for(int j = 0; j<M+1; j++)
	{
		for(int i = 0; i<N+2; i++)
		{
			v0[j][i] = 0; // Vertical velocity at n
			Rv0[j][i] = 0; // R (vertical) at n-1
		}
	}
	
	// Calculation of the constant coefficients that are used to determine the pressure
	matrix ae, aw, an, as, ap, bp;
	constant_coefficients(N+2, M+2, x, y, Sv, Sh, ae, aw, an, as, ap);
	
	// Time step (CFL condition)
	double resta = 1;
	double dtd = 0.2*rho*pow(x[2]-xvc[1],2)/mu;
	double dtc = 0.35*fabs(x[2]-xvc[1])/uref;
	double dt = min(dtd, dtc);
	
	staggx up; // Intermediate velocities
	staggy vp;
	
	
	while(resta>delta)
	{
		// STEP 1 !!! : INTERMEDIATE VELOCITY
		cout<<"Intermediate Velocity	";
		intermediate_velocities (N, M, rho, mu, method, delta, dt, x, y, xvc, yvc, Sh, Sv, V, u0, v0, Ru0, Rv0, Ru, Rv, up, vp);

		
		// STEP 2 !!! : PRESSURE
		cout<<"Pressure	";
		bp_coefficient (N+2, M+2, rho, dt, Sh, Sv, up, vp, bp);
		Gauss_Seidel (ap, aw, ae, as, an, bp, fr, delta, N+2, M+2, p);
		
		
		// STEP 3 !!! : VELOCITY
		cout<<"Velocity	";
		velocities (N, M, rho, dt, uref, x, y, p, up, vp, u, v);
		
		
		// STEP 4 !!! : TIME STEP
		cout<<"Time Step	";
		dt = time_step (dtd, x, y, u, v);
		
		
		// Comprovation
		resta = 0;
		for(int i = 0; i<N+1; i++)
		{
			for(int j = 0; j<M+2; j++)
			{
				resta = max(resta, fabs(u[j][i]-u0[j][i])/dt);
			}
		}
		for(int i = 0; i<N+2; i++)
		{
			for(int j = 0; j<M+1; j++)
			{
				resta = max(resta, fabs(v[j][i]-v0[j][i])/dt);
			}
		}
		cout<<resta<<endl;
		
		// New time step
		for(int i = 0; i<N+1; i++)
		{
			for(int j = 0; j<M+2; j++)
			{
				u0[j][i] = u[j][i];
				Ru0[j][i] = Ru[j][i];
			}
		}
		for(int i = 0; i<N+2; i++)
		{
			for(int j = 0; j<M+1; j++)
			{
				v0[j][i] = v[j][i];
				Rv0[j][i] = Rv[j][i];
			}
		}
	}
    
    output_files (x, y, u, v);
	
	return 0;
}



void coordinates(float L, int N, float xvc[], float x[])
{
	float dx = L/N;
	xvc[0] = 0;
	x[0] = 0;
	for(int i = 0; i<N; i++)
	{
		xvc[i+1] = xvc[i]+dx;
		x[i+1] = (xvc[i+1]+xvc[i])/2;
	}
	x[N+1] = L;
}


void surface(float *yvc, int M, float Sv[])
{
	for(int j = 0; j<M-1; j++)
	{
		Sv[j+1] = fabs(yvc[j]-yvc[j+1]);
	}
	Sv[0] = 0;
	Sv[M-1] = 0;
}


void volume(float *xvc, float *yvc, int N, int M, matrix& V)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(i==N-1 || j==M-1)
			{
				V[j][i] = 0;
			}
			else
			{
				V[j][i] = fabs(xvc[i]-xvc[i-1])*fabs(yvc[j]-yvc[j-1]);
			}
		}
	}
}


void constant_coefficients(int N, int M, float *x, float *y, float *Sv, float *Sh, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(j==M-1 && i!=0 && i!=N-1)
			{
				ae[j][i] = 0;
				aw[j][i] = 0;
				an[j][i] = 0;
				as[j][i] = 1;
				ap[j][i] = 1;
			}
			else if(i==0 && j==0)
			{
				ae[j][i] = 1;
				aw[j][i] = 0;
				an[j][i] = 1;
				as[j][i] = 0;
				ap[j][i] = 1;
			}
			else if(i==0 && j==M-1)
			{
				ae[j][i] = 1;
				aw[j][i] = 0;
				an[j][i] = 0;
				as[j][i] = 1;
				ap[j][i] = 1;
			}
			else if(i==0 && j!=0 && j!=M-1)
			{
				ae[j][i] = 1;
				aw[j][i] = 0;
				an[j][i] = 0;
				as[j][i] = 0;
				ap[j][i] = 1;
			}
			else if(i==N-1 && j==0)
			{
				ae[j][i] = 0;
				aw[j][i] = 1;
				an[j][i] = 1;
				as[j][i] = 0;
				ap[j][i] = 1;
			}
			else if(i==N-1 && j==M-1)
			{
				ae[j][i] = 0;
				aw[j][i] = 1;
				an[j][i] = 0;
				as[j][i] = 1;
				ap[j][i] = 1;
			}
			else if(i==N-1 && j!=0 && j!=M-1)
			{
				ae[j][i] = 0;
				aw[j][i] = 1;
				an[j][i] = 0;
				as[j][i] = 0;
				ap[j][i] = 1;
			}
			else if(j==0 && i!=0 && i!=N-1)
			{
				ae[j][i] = 0;
				aw[j][i] = 0;
				an[j][i] = 1;
				as[j][i] = 0;
				ap[j][i] = 1;
			}
			else
			{
				ae[j][i] = Sv[j]/fabs(x[i+1]-x[i]);
				aw[j][i] = Sv[j]/fabs(x[i]-x[i-1]);
				an[j][i] = Sh[i]/fabs(y[j+1]-y[j]);
				as[j][i] = Sh[i]/fabs(y[j]-y[j-1]);
				ap[j][i] = ae[j][i]+aw[j][i]+an[j][i]+as[j][i];
			}
		}
	}
}


double convective_term (string method, float xf, float x2, float x3, double u2, double u3)
{
	// 2 refers to node P, 3 to node E
	double u;
	
	if(method=="CDS")
	{
		u = u2+fabs(x2-xf)*(u3-u2)/fabs(x3-x2);
	}
	else if(method=="UDS")
	{
		int fe;
		if(u3>0)
		{
			fe = 0;
		}
		else
		{
			fe = 1;
		}
		u = u2+fe*(u3-u2);
	}
	return u;
}


void intermediate_velocities (int N, int M, float rho, float mu, string method, float delta, float dt, float* x, float* y, float *xvc, float* yvc, float* Sh, float* Sv, matrix V, staggx u0, staggy v0, staggx Ru0, staggy Rv0, staggx &Ru, staggy &Rv, staggx &up, staggy &vp)
{
	double mflowe, mfloww, mflown, mflows;
	double ue, uw, un, us;
	
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+2; j++)
		{
			// Mass flow terms (rho*v*S)
			mflowe = (rho*u0[j][i+1]+rho*u0[j][i])*Sv[j]/2;
			mfloww = (rho*u0[j][i-1]+rho*u0[j][i])*Sv[j]/2;
			mflown = (rho*v0[j][i]+rho*v0[j][i+1])*Sh[i]/2;
			mflows = (rho*v0[j-1][i]+rho*v0[j-1][i+1])*Sh[i]/2;
			
			
			// HORIZONTAL
			ue = convective_term (method, x[i+1], xvc[i], xvc[i+1], u0[j][i], u0[j][i+1]);
			uw = convective_term (method, x[i], xvc[i], xvc[i-1], u0[j][i], u0[j][i-1]);
			un = convective_term (method, yvc[j], y[j], y[j+1], u0[j][i], u0[j+1][i]);
			us = convective_term (method, yvc[j-1], y[j], y[j-1], u0[j][i], u0[j-1][i]);
			
			
			// R (horizontal)
//			if(i==0 && j==0)
//			{
//				Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-(mflowe*ue+mflown*un);
//			}
//			else if(i==0 && j==M+1)
//			{
//				Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue-mflows*us);
//			}
//			else if(i==0 && j!=0 && j!=M+1)
//			{
//				Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue+mflown*un-mflows*us);
//			}
//			else if(i==N && j==0)
//			{
//				Ru[j][i] = mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-(mflown*un-mfloww*uw);
//			}
//			else if(i==N && j==M+1)
//			{
//				Ru[j][i] = -mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(-mfloww*uw-mflows*us);
//			}
//			else if(i==N && j!=0 && j!=M+1)
//			{
//				Ru[j][i] = mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflown*un-mfloww*uw-mflows*us);
//			}
//			else if(j==0 && i!=0 && i!=N)
//			{
//				Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-(mflowe*ue+mflown*un-mfloww*uw);
//			}
//			else if(j==M+1 && i!=0 && i!=N)
//			{
//				Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue-mfloww*uw-mflows*us);
//			}
			if(i==0 || i==N || j==0 || j==M+1)
			{
				Ru[j][i] = 0;
			}
			else
			{
				Ru[j][i] = (mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue+mflown*un-mfloww*uw-mflows*us))/V[j][i];
			}
			
			// Intermediate velocity (horizontal)
			up[j][i] = u0[j][i]+dt*(1.5*Ru[j][i]-0.5*Ru0[j][i])/rho;
		}
	}
	
	double ve, vw, vn, vs;
	
	for(int i = 0; i<N+2; i++)
	{
		for(int j = 0; j<M+1; j++)
		{
			// Mass flow terms (rho*v*S)
			mflowe = (rho*u0[j+1][i]+rho*u0[j][i])*Sv[j]/2;
			mfloww = (rho*u0[j+1][i-1]+rho*u0[j][i-1])*Sv[j]/2;
			mflown = (rho*v0[j][i]+rho*v0[j+1][i])*Sh[i]/2;
			mflows = (rho*v0[j][i]+rho*v0[j-1][i])*Sh[i]/2;
			
			
			// VERTICAL
			ve = convective_term (method, xvc[i], x[i], x[i+1], v0[j][i], v0[j][i+1]);
			vw = convective_term (method, xvc[i-1], x[i], x[i-1], v0[j][i], v0[j][i-1]);
			vn = convective_term (method, y[j+1], yvc[j], yvc[j+1], v0[j][i], v0[j+1][i]);
			vs = convective_term (method, y[j], yvc[j], yvc[j-1], v0[j][i], v0[j-1][i]);
			
			// R (vertical)
//			if(i==0 && j==0)
//			{
//				Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-(mflowe*ve+mflown*vn);
//			}
//			else if(i==0 && j==M)
//			{
//				Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve-mflows*vs);
//			}
//			else if(i==0 && j!=0 && j!=M)
//			{
//				Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve+mflown*vn-mflows*vs);
//			}
//			else if(i==N+1 && j==0)
//			{
//				Rv[j][i] = mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-(mflown*vn-mfloww*vw);
//			}
//			else if(i==N+1 && j==M)
//			{
//				Rv[j][i] = -mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(-mfloww*vw-mflows*vs);
//			}
//			else if(i==N+1 && j!=0 && j!=M)
//			{
//				Rv[j][i] = mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflown*vn-mfloww*vw-mflows*vs);
//			}
//			else if(j==0 && i!=0 && i!=N+1)
//			{
//				Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-(mflowe*ve+mflown*vn-mfloww*vw);
//			}
//			else if(j==M && i!=0 && i!=N+1)
//			{
//				Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve-mfloww*vw-mflows*vs);
//			}
			if(i==0 || i==N+1 || j==0 || j==M)
			{
				Rv[j][i] = 0;
			}
			else
			{
				Rv[j][i] = (mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve+mflown*vn-mfloww*vw-mflows*vs))/V[j][i];
			}
			
			// Intermediate velocity (vertical)
			vp[j][i] = v0[j][i]+dt*(1.5*Rv[j][i]-0.5*Rv0[j][i])/rho;
		}
	}
}


void bp_coefficient (int N, int M, float rho, float dt, float* Sh, float* Sv, staggx up, staggy vp, matrix bp)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(i==0 || j==0 || i==N-1 || j==M-1)
			{
				bp[j][i] = 0;
			}
			else
			{
				bp[j][i] = -(rho*up[j][i]*Sv[j]+rho*vp[j][i]*Sh[i]-rho*up[j][i-1]*Sv[j]-rho*vp[j-1][i]*Sh[i])/dt;
			}
		}
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
				if(i==0 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==0 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+an[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==0 && j!=0 && j!=M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*T[j-1][i]+an[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N-1 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+as[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N-1 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+an[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N-1 && j!=0 && j!=M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+as[j][i]*T[j-1][i]+an[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i!=0 && i!=N-1 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*T[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i!=0 && i!=N-1 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+an[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*T[j-1][i]+an[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
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


void velocities (int N, int M, float rho, float dt, float uref, float* x, float* y, matrix p, staggx up, staggy vp, staggx &u, staggy &v)
{
	// Horizontal velocity at n+1
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+2; j++)
		{
			if(i==0 || i==N || j==0)
			{
				u[j][i] = 0;
			}
			else if(j==M+1)
			{
				u[j][i] = uref;
			}
			else
			{
				u[j][i] = up[j][i]-dt*(p[j][i+1]-p[j][i])/(rho*fabs(x[i+1]-x[i]));
			}
		}
	}
	
	// Vertical velocity at n+1
	for(int i = 0; i<N+2; i++)
	{
		for(int j = 0; j<M+1; j++)
		{
			if(j==0 || j==M || i==0 || i==N+1)
			{
				v[j][i] = 0;
			}
			else
			{
				v[j][i] = vp[j][i]-dt*(p[j+1][i]-p[j][i])/(rho*fabs(y[j+1]-y[j]));
			}
		}
	}
}


double min(double a, double b)
{
	if(a>b)
	{
		return b;
	}
	else
	{
		return a;
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


double time_step (float dtd, float* x, float* y, staggx u, staggy v)
{
	float dt;
	float dtc = 100;
	
	for(int i = 1; i<N; i++)
	{
		for(int j = 1; j<M+1; j++)
		{
			dtc = min(dtc, 0.35*fabs(x[i+1]-x[i])/fabs(u[j][i]));
		}
	}
	for(int i = 1; i<N+1; i++)
	{
		for(int j = 1; j<M; j++)
		{
			dtc = min(dtc, 0.35*fabs(y[j+1]-y[j])/fabs(v[j][i]));
		}
	}
	dt = min(dtc, dtd);
	return dt;
}


void output_files (float* x, float* y, staggx u, staggy v)
{
	ofstream resultats;
    resultats.open("Resultats.dat");
	for(int j = M+1; j>=0; j--)
	{
		for(int i = 0; i<N+2; i++)
		{
			resultats<<x[i]<<"	"<<y[j]<<"	"<<u[j][i]<<endl;
		}
		resultats<<endl;
	}
	resultats.close();
	
	ofstream resvltats;
    resvltats.open("Resvltats.dat");
	for(int j = M+1; j>=0; j--)
	{
		for(int i = 0; i<N+2; i++)
		{
			resvltats<<x[i]<<"	"<<y[j]<<"	"<<v[j][i]<<endl;
		}
		resvltats<<endl;
	}
	resvltats.close();
	
	ofstream resultsu;
    resultsu.open("u.dat");
    for(int i = M+1; i>=0; i--)
    {
    	resultsu<<y[i]<<"	"<<u[i][N/2]<<endl;
	}
    resultsu.close();
	
	ofstream resultsv;
    resultsv.open("v.dat");
    for(int i = N+1; i>=0; i--)
    {
    	resultsv<<x[i]<<"	"<<v[M/2][i]<<endl;
	}
    resultsv.close();
}
