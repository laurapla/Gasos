#include<iostream>
#include<math.h>

using namespace std;

// Numerical parameters
const int N = 10;
const int M = 10;

typedef double matrix[M][N];
typedef double staggx[M][N+1];
typedef double staggy[M+1][N];

void coordinates(float dx, int N, float xvc[], float x[]);
void surface(float *yvc, int M, float Sv[]);
void volume(float *xvc, float *yvc, int N, int M, matrix& V);
void constant_coefficients(int N, int M, float *x, float *y, float *Sv, float *Sh, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap);
double convective_term (string method, float delta, float xf, float x1, float x2, float x3, float x4, double u1, double u2, double u3, double u4);
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float* x, float fr, float delta, int N, int M, matrix& T);
double min(double a, double b);
double max(double a, double b);
void output_matrix(int N, int M, matrix mat);


int main()
{
	int Re = 100; // Reynolds number
	float L = 1; // Length of the cavity
	float rho = 1; // Density
	float uref = 1; // Reference velocity
	float mu = rho*uref*L/Re; // Viscosity
	
	
	string method = "CDS";
	float delta = 0.001; // Precision of the simulation
	float fr = 1; // Relaxation factor
	
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
	matrix p; // Values in the nodes
	staggx u, u0, u00, uant, Ru, Ru0; // Values in the points given by the staggered meshes
	staggy v, v0, v00, vant, Rv, Rv0;
	
	// Inicialization
	for(int j = 0; j<M; j++)
	{
		for(int i = 0; i<N+1; i++)
		{
			if(j==M-1)
			{
				u[j][i] = uref; // Horizontal velocity
			}
			else
			{
				u[j][i] = 0; // Horizontal velocity at n+1
			}
			u0[j][i] = 0; // Horizontal velocity at n
			u00[j][i] = 0; // Horizontal velocity at n-1
			uant[j][i] = u[j][i]+10;
			Ru0[j][i] = 0;
		}
	}
	for(int j = 0; j<M+1; j++)
	{
		for(int i = 0; i<N; i++)
		{
			v[j][i] = 0; // Vertical velocity at n+1
			v0[j][i] = 0; // Vertical velocity at n
			v00[j][i] = 0; // Vertical velocity at n-1
			vant[j][i] = v[j][i]+10;
			Rv0[j][i] = 0;
		}
	}
	
	// Calculation of the constant coefficients that are used to determine the pressure
	matrix ae, aw, an, as, ap, bp;
	constant_coefficients(N, M, x, y, Sv, Sh, ae, aw, an, as, ap);
	
	double resta = 1;
	double dtd = 0.2*rho*pow(x[1]-xvc[1],2)/mu;
	double dtc = 0.35*fabs(x[1]-xvc[1])/uref;
	double dt = min(dtd, dtc);
	
	double mflowe, mfloww, mflown, mflows;
	double ue, uw, un, us; // Horizontal velocities
	double ve, vw, vn, vs; // Vertical velocities
	matrix up, vp; // Intermediate velocities
	
	
	while(resta>delta)
	{
		// STEP 1 !!! : INTERMEDIATE VELOCITY
		cout<<"Intermediate Velocity	";
		for(int i = 0; i<N+1; i++)
		{
			for(int j = 0; j<M; j++)
			{
				// Mass flow terms (rho*v*S)
				mflowe = (rho*u0[j][i+1]+rho*u0[j][i])*Sv[j]/2;
				mfloww = (rho*u0[j][i-1]+rho*u0[j][i])*Sv[j]/2;
				mflown = (rho*v0[j][i]+rho*v0[j+1][i])*Sh[i]/2;
				mflows = (rho*v0[j][i]+rho*v0[j-1][i])*Sh[i]/2;
				
				
				// HORIZONTAL
				ue = convective_term (method, delta, x[i], xvc[i-1], xvc[i], xvc[i+1], xvc[i+2], u0[j][i-1], u0[j][i], u0[j][i+1], u0[j][i+2]);
				uw = convective_term (method, delta, x[i-1], xvc[i+1], xvc[i], xvc[i-1], xvc[i-2], u0[j][i+1], u0[j][i], u0[j][i-1], u0[j][i-2]);
				un = convective_term (method, delta, y[j], yvc[j-1], yvc[j], yvc[j+1], yvc[j+2], u0[j-1][i], u0[j][i], u0[j+1][i], u0[j+2][i]);
				us = convective_term (method, delta, y[j-1], yvc[j+1], yvc[j], yvc[j-1], yvc[j-2], u0[j+1][i], u0[j][i], u0[j-1][i], u0[j-2][i]);
				
				
				// R (horizontal)
				if(i==0 && j==0)
				{
					Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-(mflowe*ue+mflown*un);
				}
				else if(i==0 && j==M-1)
				{
					Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue-mflows*us);
				}
				else if(i==0 && j!=0 && j!=M-1)
				{
					Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue+mflown*un-mflows*us);
				}
				else if(i==N && j==0)
				{
					Ru[j][i] = mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-(mflown*un-mfloww*uw);
				}
				else if(i==N && j==M-1)
				{
					Ru[j][i] = -mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(-mfloww*uw-mflows*us);
				}
				else if(i==N && j!=0 && j!=M-1)
				{
					Ru[j][i] = mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflown*un-mfloww*uw-mflows*us);
				}
				else if(j==0 && i!=0 && i!=N)
				{
					Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-(mflowe*ue+mflown*un-mfloww*uw);
				}
				else if(j==M-1 && i!=0 && i!=N)
				{
					Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue-mfloww*uw-mflows*us);
				}
				else
				{
					Ru[j][i] = mu*(u0[j][i+1]-u0[j][i])*Sv[j]/fabs(xvc[i+1]-xvc[i])+mu*(u0[j+1][i]-u0[j][i])*Sh[i]/fabs(y[j+1]-y[j])-mu*(u0[j][i]-u0[j][i-1])*Sv[j]/fabs(xvc[i]-xvc[i-1])-mu*(u0[j][i]-u0[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue+mflown*un-mfloww*uw-mflows*us);
				}
				
				// Intermediate velocity (horizontal)
				up[j][i] = u0[j][i]+dt*(1.5*Ru[j][i]-0.5*Ru0[j][i])/rho;
			}
		}
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M+1; j++)
			{
				// Mass flow terms (rho*v*S)
				mflowe = (rho*u0[j][i+1]+rho*u0[j][i])*Sv[j]/2;
				mfloww = (rho*u0[j][i-1]+rho*u0[j][i])*Sv[j]/2;
				mflown = (rho*v0[j][i]+rho*v0[j+1][i])*Sh[i]/2;
				mflows = (rho*v0[j][i]+rho*v0[j-1][i])*Sh[i]/2;
				
				
				// VERTICAL
				ve = convective_term (method, delta, xvc[i+1], x[i-1], x[i], x[i+1], x[i+2], v0[j][i-1], v0[j][i], v0[j][i+1], v0[j][i+2]);
				vw = convective_term (method, delta, xvc[i], x[i+1], x[i], x[i-1], x[i-2], v0[j][i+1], v0[j][i], v0[j][i-1], v0[j][i-2]);
				vn = convective_term (method, delta, yvc[j+1], y[j-1], y[j], y[j+1], y[j+2], v0[j-1][i], v0[j][i], v0[j+1][i], v0[j+2][i]);
				vs = convective_term (method, delta, yvc[j], y[j+1], y[j], y[j-1], y[j-2], v0[j+1][i], v0[j][i], v0[j-1][i], v0[j-2][i]);
				
				// R (vertical)
				if(i==0 && j==0)
				{
					Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-(mflowe*ve+mflown*vn);
				}
				else if(i==0 && j==M)
				{
					Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve-mflows*vs);
				}
				else if(i==0 && j!=0 && j!=M)
				{
					Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve+mflown*vn-mflows*vs);
				}
				else if(i==N-1 && j==0)
				{
					Rv[j][i] = mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-(mflown*vn-mfloww*vw);
				}
				else if(i==N-1 && j==M)
				{
					Rv[j][i] = -mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(-mfloww*vw-mflows*vs);
				}
				else if(i==N-1 && j!=0 && j!=M)
				{
					Rv[j][i] = mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflown*vn-mfloww*vw-mflows*vs);
				}
				else if(j==0 && i!=0 && i!=N-1)
				{
					Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-(mflowe*ve+mflown*vn-mfloww*vw);
				}
				else if(j==M && i!=0 && i!=N-1)
				{
					Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve-mfloww*vw-mflows*vs);
				}
				else
				{
					Rv[j][i] = mu*(v0[j][i+1]-v0[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v0[j+1][i]-v0[j][i])*Sh[i]/fabs(yvc[j+1]-yvc[j])-mu*(v0[j][i]-v0[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v0[j][i]-v0[j-1][i])*Sh[i]/fabs(yvc[j]-yvc[j-1])-(mflowe*ve+mflown*vn-mfloww*vw-mflows*vs);
				}
				
				// Intermediate velocity (vertical)
				vp[j][i] = v0[j][i]+dt*(1.5*Rv[j][i]-0.5*Rv0[j][i])/rho;
			}
		}
		
		
		// STEP 2 !!! : PRESSURE
		cout<<"Pressure	";
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M; j++)
			{
				bp[j][i] = -(rho*up[j][i+1]*Sv[j]+rho*vp[j+1][i]*Sh[i]-rho*up[j][i]*Sv[j]-rho*vp[j][i]*Sh[i])/dt;
			}
		}
		Gauss_Seidel (ap, aw, ae, as, an, bp, x, fr, delta, N, M, p);
		
		
		// STEP 3 !!! : VELOCITY
		cout<<"Velocity	";
		// Horizontal velocity
		for(int i = 0; i<N+1; i++)
		{
			for(int j = 0; j<M; j++)
			{
				if(i==0)
				{
					u[j][i] = up[j][i]-dt*p[j][i]/(rho*fabs(x[i]-xvc[i]));
				}
				else if(j==M-1)
				{
					u[j][i] = uref;
				}
				else
				{
					u[j][i] = up[j][i]-dt*(p[j][i]-p[j][i-1])/(rho*fabs(x[i]-x[i-1]));
				}
			}
		}
		
		// Vertical velocity
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M+1; j++)
			{
				if(j==0)
				{
					v[j][i] = vp[j][i]-dt*p[j][i]/(rho*fabs(y[j]-yvc[j]));
				}
				else
				{
					v[j][i] = vp[j][i]-dt*(p[j][i]-p[j-1][i])/(rho*fabs(y[j]-y[j-1]));
				}
			}
		}
		
		
		// STEP 4 !!! : TIME STEP
		cout<<"Time Step	";
		dtc = 100;
		for(int i = 0; i<N+1; i++)
		{
			for(int j = 0; j<M; j++)
			{
				dtc = min(dtc, 0.35*fabs(x[3]-x[2])/fabs(u[j][i]));
			}
		}
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M+1; j++)
			{
				dtc = min(dtc, 0.35*fabs(y[3]-y[2])/fabs(v[j][i]));
			}
		}
		dt = min(dtc, dtd);
		
		
		// Comprovation
		resta = 0;
		for(int i = 0; i<N+1; i++)
		{
			for(int j = 0; j<M; j++)
			{
				resta = max(resta, fabs(u[j][i]-u0[j][i]));
			}
		}
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M+1; j++)
			{
				resta = max(resta, fabs(v[j][i]-v0[j][i]));
			}
		}
		cout<<resta<<endl;
		
		// New time step
		for(int i = 0; i<N+1; i++)
		{
			for(int j = 0; j<M; j++)
			{
				u0[j][i] = u[j][i];
				Ru0[j][i] = Ru[j][i];
			}
		}
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M+1; j++)
			{
				v0[j][i] = v[j][i];
				Rv0[j][i] = Rv[j][i];
			}
		}
	}
	
	
//	for(int j = M-1; j>=0; j--)
//	{
//		for(int i = 0; i<N+1; i++)
//		{
//			cout<<u[j][i]<<"	";
//		}
//		cout<<endl;
//	}
//	for(int i = 0; i<N; i++)
//	{
//		cout<<x[i]<<endl;
//	}
//	output_matrix(N, M, p);
	
	
	
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


void constant_coefficients(int N, int M, float *x, float *y, float *Sv, float *Sh, matrix& ae, matrix& aw, matrix& an, matrix& as, matrix& ap)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(j==M-1)
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


double convective_term (string method, float delta, float xf, float x1, float x2, float x3, float x4, double u1, double u2, double u3, double u4)
{
	// 1 refers to node W, 2 to node P, 3 to node E, and 4 to node EE
	double u, ud, uu, uuu;
	float xd, xu, xuu;
	float resta = 1;
	
	if(method=="CDS")
	{
		u = 0.5*(u1+u2);
	}
	else if(method=="UDS")
	{
		int fe;
		if(u2>0)
		{
			fe = 0;
		}
		else
		{
			fe = 1;
		}
		u = u1+fe*(u2-u1);
	}
	else if(method=="QUICK")
	{
		float g1, g2;
		if(u2>0)
		{
			ud = u3;
			uu = u2;
			uuu = u1;
			xd = x3;
			xu = x2;
			xuu = x1;
		}
		else
		{
			ud = u2;
			uu = u3;
			uuu = u4;
			xd = x2;
			xu = x3;
			xuu = x4;
		}
		g1 = (xf-xu)*(xf-xuu)/((xd-xu)*(xd-xuu));
		g2 = (xf-xu)*(xd-xf)/((xu-xuu)*(xd-xuu));
		u = uu+g1*(ud-uu)+g2*(uu-uuu);
	}
	return u;
}


// Solver (using Gauss-Seidel)
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float* x, float fr, float delta, int N, int M, matrix& T)
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
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==0 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==0 && j!=0 && j!=M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j-1][i]+an[j][i]*T[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N-1 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+as[j][i]*Tcalc[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N-1 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+an[j][i]*T[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i==N && j!=0 && j!=M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+as[j][i]*Tcalc[j-1][i]+an[j][i]*T[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i!=0 && i!=N-1 && j==M-1)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else if(i!=0 && i!=N-1 && j==0)
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
				}
				else
				{
					T[j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j-1][i]+an[j][i]*T[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
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


void output_matrix(int N, int M, matrix mat)
{
	for(int j = M-1; j>=0; j--)
	{
		for(int i = 0; i<N; i++)
		{
			cout<<mat[j][i]<<"	";
		}
		cout<<endl;
	}
}
