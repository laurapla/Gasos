#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

// Numerical parameters
const int N = 50;
const int M = 50;

typedef double matrix[M+2][N+2];
typedef double staggx[M+2][N+1];
typedef double staggy[M+1][N+2];

void coordinates(float L, int N, double xvc[], double x[]);
void surface(double *yvc, int M, double Sv[]);
void volume(double *xvc, double *yvc, int N, int M, matrix& V);
void initial_conditions(int N, int M, staggx& u0, staggx& Ru0, staggy& v0, staggy& Rv0, matrix& T0);
void constant_coefficients(int N, int M, double *x, double *y, double *Sv, double *Sh, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap);
void temperature_coefficients(int N, int M, double dt, double* x, double* y, double* Sv, double* Sh, matrix V, staggx u, staggy v, matrix T0, matrix &aTe, matrix &aTw, matrix &aTn, matrix &aTs, matrix &aTp, matrix &bTp);
double convective_term (double xf, double x2, double x3, double u2, double u3);
void intermediate_velocities (int N, int M, float Pr, int Ra, double dt, double* x, double* y, double *xvc, double* yvc, double* Sh, double* Sv, matrix V, matrix T0, staggx u0, staggy v0, staggx Ru0, staggy Rv0, staggx &Ru, staggy &Rv, staggx &up, staggy &vp);
void bp_coefficient (int N, int M, double dt, double* Sh, double* Sv, staggx up, staggy vp, matrix bp);
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float fr, float delta, int N, int M, matrix& T);
void velocities (int N, int M, double dt, double* x, double* y, matrix p, staggx up, staggy vp, staggx &u, staggy &v);
double min(double a, double b);
double max(double a, double b);
double time_step (double dtd, double* x, double* y, staggx u, staggy v);
double error (int N, int M, staggx u, staggy v, staggx u0, staggy v0, matrix T, matrix T0);
void heat_flux(int N, int M, double* x, staggx u, matrix T, matrix Q);
void Nusselt(int N, int M, double* x, double* yvc, matrix Q, double Nu[]);
void maximum_planes (int N, int M, double* x, double* y, staggx u, staggy v);
void output_files (int N, int M, float L, double* x, double* y, double* xvc, double* yvc, staggx u, staggy v, matrix T, double* Nu);


int main()
{
	float Pr = 0.71; // Prandtl number
	int Ra = 1e6; // Rayleigh number
	float L = 1; // Length of the cavity
	
	float delta = 1e-4; // Precision of the simulation
	float fr = 1.2; // Relaxation factor
	
	cout<<"Program started"<<endl;
	cout<<"Pr="<<Pr<<endl<<endl;
	cout<<"Ra="<<Ra<<endl<<endl;
	
	// Coordinates
	double xvc[N+1], yvc[M+1], x[N+2], y[M+2];
	coordinates(L, N, xvc, x);
	coordinates(L, M, yvc, y);
	
	// Surfaces
	double Sh[N+2], Sv[M+2];
	matrix V;
	surface(xvc, N+2, Sh); // Horizontal surface
	surface(yvc, M+2, Sv); // Vertical surface
	volume(xvc, yvc, N+2, M+2, V); // Volume
	
	
	// Properties that are going to be calculated
	matrix p, T, T0, Q; // Values in the nodes (pressure)
	double Nu[N+2]; // Nusselt number
	staggx u, u0, Ru0; // Values in the points given by the staggered meshes (velocities)
	staggy v, v0, Rv0;
	
	// Inicialization
	initial_conditions(N, M, u0, Ru0, v0, Rv0, T0);
	
	matrix aTe, aTw, aTn, aTs, aTp, bTp;
	
	// Calculation of the constant coefficients that are used to determine the pressure
	matrix ae, aw, an, as, ap, bp;
	constant_coefficients(N+2, M+2, x, y, Sv, Sh, ae, aw, an, as, ap);
	
	// Time step (CFL condition)
	double resta = 1;
	double dtd = 0.2*pow(x[2]-xvc[1],2)/Pr;
	double dtc = 0.35*fabs(x[2]-xvc[1]);
	double dt = min(dtd, dtc);
	
	staggx up, Ru; // Intermediate velocities
	staggy vp, Rv;
	
	cout<<"Solving..."<<endl;
	// Fractional Step Method
	while(resta>delta)
	{
		// STEP 1: INTERMEDIATE VELOCITY
		intermediate_velocities (N, M, Pr, Ra, dt, x, y, xvc, yvc, Sh, Sv, V, T0, u0, v0, Ru0, Rv0, Ru, Rv, up, vp);

		
		// STEP 2: PRESSURE
		bp_coefficient (N+2, M+2, dt, Sh, Sv, up, vp, bp);
		Gauss_Seidel (ap, aw, ae, as, an, bp, fr, delta, N+2, M+2, p);
		
		
		// STEP 3: VELOCITY
		velocities (N, M, dt, x, y, p, up, vp, u, v);
		
		
		// STEP 4: TEMPERATURE
		temperature_coefficients(N, M, dt, x, y, Sv, Sh, V, u, v, T0, aTe, aTw, aTn, aTs, aTp, bTp);
		Gauss_Seidel (aTp, aTw, aTe, aTs, aTn, bTp, fr, delta, N+2, M+2, T);
		
		
		// STEP 5: TIME STEP
		dt = time_step (dtd, x, y, u, v);
		
		
		// Comprovation
		resta = error (N, M, u, v, u0, v0, T, T0);
		
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
		for(int j = 0; j<M+2; j++)
		{
			for(int i = 0; i<N+2; i++)
			{
				T0[j][i] = T[j][i];
			}
		}
	}
	
	// Results
	heat_flux(N, M, x, u, T, Q);
	Nusselt(N, M, x, y, Q, Nu);
    cout<<endl<<"Creating some output files..."<<endl;
    output_files (N, M, L, x, y, xvc, yvc, u, v, T, Nu);
	maximum_planes (N, M, x, y, u, v);
	
	return 0;
}


// Coordinates of the control volumes (x -> nodes, xvc -> faces)
void coordinates(float L, int N, double xvc[], double x[])
{
	double dx = L/N;
	xvc[0] = 0;
	x[0] = 0;
	for(int i = 0; i<N; i++)
	{
		xvc[i+1] = xvc[i]+dx;
		x[i+1] = (xvc[i+1]+xvc[i])/2;
	}
	x[N+1] = L;
}


// Surfaces of the control volumes
void surface(double *yvc, int M, double Sv[])
{
	for(int j = 0; j<M-1; j++)
	{
		Sv[j+1] = fabs(yvc[j]-yvc[j+1]);
	}
	Sv[0] = 0;
	Sv[M-1] = 0;
}


// Volume of each control volume
void volume(double *xvc, double *yvc, int N, int M, matrix& V)
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


// Initial conditions
void initial_conditions(int N, int M, staggx& u0, staggx& Ru0, staggy& v0, staggy& Rv0, matrix& T0)
{
	for(int j = 0; j<M+2; j++)
	{
		for(int i = 0; i<N+1; i++)
		{
			u0[j][i] = 0; // Horizontal velocity at n
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
	for(int j = 0; j<M+2; j++)
	{
		for(int i = 0; i<N+2; i++)
		{
			if(i==0)
			{
				T0[j][i] = 1;
			}
			else
			{
				T0[j][i] = 0;
			}
		}
	}
}


// Calculation of the constant coefficients (ae, aw, an, as, ap) of the Poisson equation (pressure)
void constant_coefficients(int N, int M, double *x, double *y, double *Sv, double *Sh, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap)
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


// Coefficients used to calculate the temperature
void temperature_coefficients(int N, int M, double dt, double* x, double* y, double* Sv, double* Sh, matrix V, staggx u, staggy v, matrix T0, matrix &aTe, matrix &aTw, matrix &aTn, matrix &aTs, matrix &aTp, matrix &bTp)
{
	double Fe, Fw, Fn, Fs;
	double De, Dw, Dn, Ds;
	
	for(int j = 0; j<M+2; j++)
	{
		for(int i = 0; i<N+2; i++)
		{
			if(i==0)
			{
				aTe[j][i] = 0;
				aTw[j][i] = 0;
				aTn[j][i] = 0;
				aTs[j][i] = 0;
				aTp[j][i] = 1;
				bTp[j][i] = 1;
			}
			else if(i==N+1)
			{
				aTe[j][i] = 0;
				aTw[j][i] = 0;
				aTn[j][i] = 0;
				aTs[j][i] = 0;
				aTp[j][i] = 1;
				bTp[j][i] = 0;
			}
			else if(j==0 && i!=0 && i!=N+1)
			{
				aTe[j][i] = 0;
				aTw[j][i] = 0;
				aTn[j][i] = 1;
				aTs[j][i] = 0;
				aTp[j][i] = 1;
				bTp[j][i] = 0;
			}
			else if(j==M+1 && i!=0 && i!=N+1)
			{
				aTe[j][i] = 0;
				aTw[j][i] = 0;
				aTn[j][i] = 0;
				aTs[j][i] = 1;
				aTp[j][i] = 1;
				bTp[j][i] = 0;
			}
			else
			{
				// Mass flow terms (v*S)
				Fe = u[j][i]*Sv[j];
				Fw = u[j][i-1]*Sv[j];
				Fn = v[j][i]*Sh[i];
				Fs = v[j-1][i]*Sh[i];
				
				// Areas and distances
				De = Sv[j]/fabs(x[i+1]-x[i]);
				Dw = Sv[j]/fabs(x[i]-x[i-1]);
				Dn = Sh[i]/fabs(y[j+1]-y[j]);
				Ds = Sh[i]/fabs(y[j+1]-y[j]);
			
				aTe[j][i] = De-0.5*Fe;
				aTw[j][i] = Dw+0.5*Fw;
				aTn[j][i] = Dn-0.5*Fn;
				aTs[j][i] = Ds+0.5*Fs;
				aTp[j][i] = aTe[j][i]+aTw[j][i]+aTn[j][i]+aTs[j][i]+V[j][i]/dt;
				bTp[j][i] = T0[j][i]*V[j][i]/dt;
			}
		}
	}
}


// Computation of the velocity in the convective term using CDS
double convective_term (double xf, double x2, double x3, double u2, double u3)
{
	// 2 refers to node P, 3 to node E
	double u;
	u = u2+fabs(x2-xf)*(u3-u2)/fabs(x3-x2);
	
	return u;
}


// Calculation of the intermediate velocities
void intermediate_velocities (int N, int M, float Pr, int Ra, double dt, double* x, double* y, double *xvc, double* yvc, double* Sh, double* Sv, matrix V, matrix T0, staggx u0, staggy v0, staggx Ru0, staggy Rv0, staggx &Ru, staggy &Rv, staggx &up, staggy &vp)
{
	double mflowe, mfloww, mflown, mflows;
	double ue, uw, un, us;
	
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+2; j++)
		{
			// Mass flow terms (v*S)
			mflowe = (u0[j][i+1]+u0[j][i])*Sv[j]/2;
			mfloww = (u0[j][i-1]+u0[j][i])*Sv[j]/2;
			mflown = (v0[j][i]+v0[j][i+1])*Sh[i]/2;
			mflows = (v0[j-1][i]+v0[j-1][i+1])*Sh[i]/2;
			
			
			// HORIZONTAL
			ue = convective_term (x[i+1], xvc[i], xvc[i+1], u0[j][i], u0[j][i+1]);
			uw = convective_term (x[i], xvc[i], xvc[i-1], u0[j][i], u0[j][i-1]);
			un = convective_term (yvc[j], y[j], y[j+1], u0[j][i], u0[j+1][i]);
			us = convective_term (yvc[j-1], y[j], y[j-1], u0[j][i], u0[j-1][i]);
			
			
			// R (horizontal)
			if(i==0 || i==N || j==0 || j==M+1)
			{
				Ru[j][i] = 0;
			}
			else
			{
				Ru[j][i] = (Pr*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+Pr*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-Pr*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-Pr*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue+mflown*un-mfloww*uw-mflows*us))/V[j][i];
			}
			
			// Intermediate velocity (horizontal)
			up[j][i] = u0[j][i]+dt*(1.5*Ru[j][i]-0.5*Ru0[j][i]);
		}
	}
	
	double ve, vw, vn, vs;
	
	for(int i = 0; i<N+2; i++)
	{
		for(int j = 0; j<M+1; j++)
		{
			// Mass flow terms (v*S)
			mflowe = (u0[j+1][i]+u0[j][i])*Sv[j]/2;
			mfloww = (u0[j+1][i-1]+u0[j][i-1])*Sv[j]/2;
			mflown = (v0[j][i]+v0[j+1][i])*Sh[i]/2;
			mflows = (v0[j][i]+v0[j-1][i])*Sh[i]/2;
			
			
			// VERTICAL
			ve = convective_term (xvc[i], x[i], x[i+1], v0[j][i], v0[j][i+1]);
			vw = convective_term (xvc[i-1], x[i], x[i-1], v0[j][i], v0[j][i-1]);
			vn = convective_term (y[j+1], yvc[j], yvc[j+1], v0[j][i], v0[j+1][i]);
			vs = convective_term (y[j], yvc[j], yvc[j-1], v0[j][i], v0[j-1][i]);
			
			// R (vertical)
			if(i==0 || i==N+1 || j==0 || j==M)
			{
				Rv[j][i] = 0;
			}
			else
			{
				Rv[j][i] = (Pr*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+Pr*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-Pr*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-Pr*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve+mflown*vn-mfloww*vw-mflows*vs))/V[j][i]+Pr*Ra*(T0[j][i]+T0[j+1][i])/2;
			}
			
			// Intermediate velocity (vertical)
			vp[j][i] = v0[j][i]+dt*(1.5*Rv[j][i]-0.5*Rv0[j][i]);
		}
	}
}


// Calculation of the bp coefficient of the Poisson equation (pressure)
void bp_coefficient (int N, int M, double dt, double* Sh, double* Sv, staggx up, staggy vp, matrix bp)
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
				bp[j][i] = -(up[j][i]*Sv[j]+vp[j][i]*Sh[i]-up[j][i-1]*Sv[j]-vp[j-1][i]*Sh[i])/dt;
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


// Calculation of the velocity with the pressure correction
void velocities (int N, int M, double dt, double* x, double* y, matrix p, staggx up, staggy vp, staggx &u, staggy &v)
{
	// Horizontal velocity at n+1
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+2; j++)
		{
			if(i==0 || i==N || j==0 || j==M+1)
			{
				u[j][i] = 0;
			}
			else
			{
				u[j][i] = up[j][i]-dt*(p[j][i+1]-p[j][i])/(fabs(x[i+1]-x[i]));
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
				v[j][i] = vp[j][i]-dt*(p[j+1][i]-p[j][i])/(fabs(y[j+1]-y[j]));
			}
		}
	}
}


// Returns the minimum value
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


// Returns the maximum value
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


// Calculation of the proper time step (CFL condition)
double time_step (double dtd, double* x, double* y, staggx u, staggy v)
{
	double dt;
	double dtc = 100;
	
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


// Difference between the previous and the actual time step
double error (int N, int M, staggx u, staggy v, staggx u0, staggy v0, matrix T, matrix T0)
{
	double resta = 0;
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+2; j++)
		{
			resta = max(resta, fabs(u[j][i]-u0[j][i]));
		}
	}
	for(int i = 0; i<N+2; i++)
	{
		for(int j = 0; j<M+1; j++)
		{
			resta = max(resta, fabs(v[j][i]-v0[j][i]));
		}
	}
	for(int j = 0; j<M+2; j++)
	{
		for (int i = 0; i<N+2; i++)
		{
			resta = max(resta, fabs(T[j][i]-T0[j][i]));
		}
	}
	
	cout<<resta<<endl;
	return resta;
}


// Heat flux in the horizontal direction at any point in the cavity
void heat_flux(int N, int M, double* x, staggx u, matrix T, matrix Q)
{
	for(int j = 0; j<M+2; j++)
	{
		for(int i = 0; i<N+1; i++)
		{
			if(i==0)
			{
				Q[j][i] = u[j][i]*T[j][i]-(T[j][i+1]-T[j][i])/fabs(x[i+1]-x[i]);
			}
			else if(i==N)
			{
				Q[j][i] = Q[j][0];
			}
			else
			{
				Q[j][i] = 0.5*u[j][i]*(T[j][i]+T[j][i+1])-(T[j][i+1]-T[j][i])/fabs(x[i+1]-x[i]);
			}
		}
	}
}


// Computation of the Nusselt numbers
void Nusselt(int N, int M, double* x, double* yvc, matrix Q, double Nu[])
{
	for(int i = 0; i<N+1; i++)
	{
		Nu[i] = 0;
		for(int j = 0; j<M+1; j++)
		{
			Nu[i] = Nu[i]+(yvc[j+1]-yvc[j])*Q[j][i];
		}
		
	}
	
	double Numax = -100;
	double Numin = 100;
	double Nuavg = 0;
	double Nu0 = Nu[0];
	double Nu12 = (Nu[N/2+1]+Nu[N/2+2])/2;
	int jmax, jmin;
	
	for(int i = 0; i<N+2; i++)
	{
		Nuavg = Nuavg+(x[i+1]-x[i])*Nu[i];
	}
	for(int j = 0; j<M+2; j++)
	{
		if(Q[j][0]>Numax)
		{
			Numax = Q[j][0];
			jmax = j;
		}
		if(Q[j][0]<Numin)
		{
			Numin = Q[j][0];
			jmin = j;
		}
	}
	cout<<endl<<endl;
	cout<<"Nu average = "<<Nuavg<<endl;
	cout<<"Nu0 = "<<Nu0<<endl;
	cout<<"Nu1/2 = "<<Nu12<<endl;
	cout<<"Nu max = "<<Numax<<" at y = "<<yvc[jmax]<<endl;
	cout<<"Nu min = "<<Numin<<" at y = "<<yvc[jmin]<<endl;
}


// Maximum velocity at the central horizontal and vertical planes
void maximum_planes (int N, int M, double* x, double* y, staggx u, staggy v)
{
	double umax = 0, vmax = 0;
	int imax, jmax;
	double uavg, vavg;
	for(int j = 0; j<M+2; j++)
	{
		uavg = (u[j][N/2+1]+u[j][N/2])/2;
		if(uavg>umax)
		{
			umax = uavg;
			jmax = j;
		}
	}
	for(int i = 0; i<N+2; i++)
	{
		vavg = (v[M/2+1][i]+v[M/2][i])/2;
		if(vavg>vmax)
		{
			vmax = vavg;
			imax = i;
		}
	}
	
	cout<<"u max = "<<umax<<" at y = "<<y[jmax]<<endl;
	cout<<"v max = "<<vmax<<" at x = "<<x[imax]<<endl;
}


// Output of the results
void output_files (int N, int M, float L, double* x, double* y, double* xvc, double* yvc, staggx u, staggy v, matrix T, double* Nu)
{
	// Horizontal coordinates
    ofstream xx;
	xx.open("x.dat");
    for(int i = 0; i<N+2; i++)
    {
    	xx<<x[i]<<endl;
	}
    xx.close();
    
    // Vertical coordinates
	ofstream yy;
    yy.open("y.dat");
    for(int j = M+1; j>=0; j--)
    {
    	yy<<y[j]<<endl;
	}
    yy.close();
	
	// Horizontal velocities
	ofstream resultats;
    resultats.open("Resultats.dat");
	for(int i = 0; i<N+1; i++)
	{
		for(int j = 0; j<M+2; j++)
		{
			resultats<<xvc[i]<<"	"<<y[j]<<"	"<<u[j][i]<<endl;
		}
		resultats<<endl;
	}
	resultats.close();
	
	// Vertical velocities
	ofstream resvltats;
    resvltats.open("Resvltats.dat");
	for(int i = 0; i<N+2; i++)
	{
		for(int j = 0; j<M+1; j++)
		{
			resvltats<<x[i]<<"	"<<yvc[j]<<"	"<<v[j][i]<<endl;
		}
		resvltats<<endl;
	}
	resvltats.close();
    
    // Matrix of horizontal velocities
	ofstream result;
    result.open("Matrixu.dat");
	for(int j = M+1; j>=0; j--)
	{
		for(int i = 0; i<N+2; i++)
		{
			if(i==0 || i==N+1)
			{
				result<<u[j][i]<<"	";
			}
			else
			{
				result<<convective_term (x[i], xvc[i-1], xvc[i], u[j][i-1], u[j][i])<<"	";
			}
		}
		result<<endl;
	}
	result.close();
	
	// Matrix of vertical velocities
	ofstream resvlt;
    resvlt.open("Matrixv.dat");
	for(int j = M+1; j>=0; j--)
	{
		for(int i = 0; i<N+2; i++)
		{
			if(j==0 && j==M+1)
			{
				resvlt<<v[j][i]<<"	";
			}
			else
			{
				resvlt<<convective_term (y[j], yvc[j-1], yvc[j], v[j-1][i], v[j][i])<<"	";
			}
		}
		resvlt<<endl;
	}
	resvlt.close();
    
    
    // Temperature
    ofstream temperature;
    temperature.open("Temperatura.dat");
    for(int j = M+1; j>=0; j--)
    {
    	for(int i = 0; i<N+2; i++)
    	{
    		temperature<<T[j][i]<<"	";
		}
		temperature<<endl;
	}
	temperature.close();
	
	// Nusselt number
	ofstream nuss;
	nuss.open("Nusselt.dat");
	for(int i = 0; i<N+1; i++)
	{
		nuss<<xvc[i]<<"	"<<Nu[i]<<endl;
	}
	nuss.close();
}
