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
void output_matrix(int N, int M, matrix mat);

int main()
{
	int Re = 100; // Reynolds number
	float L = 1; // Length of the cavity
	float rho = 1; // Density
	float mu = 0.00001; // Viscosity
	float uref = Re*mu/(rho*L); // Reference velocity
	
	string method = "UDS";
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
	matrix p, u, v, u0, v0; // Values in the nodes
	staggx ustagg; // Values in the points given by the staggered meshes
	staggy vstagg;
	
	// Inicialization
	for(int j = 0; j<M; j++)
	{
		for(int i = 0; i<N; i++)
		{
			u[j][i] = 0; // Horizontal velocity
			u0[j][i] = u[j][i];
			ustagg[j][i] = u[j][i];
			v[j][i] = 0; // Vertical velocity
			v0[j][i] = v[j][i];
			vstagg[j][i] = v[j][i];
		}
	}
	
	
	// Calculation of the constant coefficients that are used to determine the pressure
	matrix ae, aw, an, as, ap, bp;
	constant_coefficients(N, M, x, y, Sv, Sh, ae, aw, an, as, ap);
	
	double resta = 1;
	float dt = 0.01;
	
	double mflowe, mfloww, mflown, mflows;
	double ue, uw, un, us; // Horizontal velocities
	double ve, vw, vn, vs; // Vertical velocities
	matrix up, vp; // Intermediate velocities
	double Ru, Ru0, Rv, Rv0;
	
	while(resta>delta)
	{
		// STEP 1 !!!
		
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M; j++)
			{
				// Mass flow terms (rho*v*S)
				mflowe = (rho*u[j][i+1]+rho*u[j][i])*Sv[j]/2;
				mfloww = (rho*u[j][i-1]+rho*u[j][i])*Sv[j]/2;
				mflown = (rho*v[j][i]+rho*v[j+1][i])*Sh[i]/2;
				mflows = (rho*v[j][i]+rho*v[j-1][i])*Sh[i]/2;
				
				
				// HORIZONTAL
				ue = convective_term (method, delta, xvc[i+1], x[i-1], x[i], x[i+1], x[i+2], u[j][i-1], u[j][i], u[j][i+1], u[j][i+2]);
				uw = convective_term (method, delta, xvc[i], x[i+1], x[i], x[i-1], x[i-2], u[j][i+1], u[j][i], u[j][i-1], u[j][i-2]);
				un = convective_term (method, delta, yvc[j+1], y[j-1], y[j], y[j+1], y[j+2], u[j-1][i], u[j][i], u[j+1][i], u[j+2][i]);
				us = convective_term (method, delta, yvc[j], y[j+1], y[j], y[j-1], y[j-2], u[j+1][i], u[j][i], u[j-1][i], u[j-2][i]);
				
				// R (horizontal)
				Ru = mu*(u[j][i+1]-u[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(u[j+1][i]-u[j][i])*Sh[j]/fabs(y[j+1]-y[j])-mu*(u[j][i]-u[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(u[j][i]-u[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ue+mflown*un-mfloww*uw-mflows*us);
				Ru0 = 0;
				
				// Intermediate velocity (horizontal)
				up[j][i] = u[j][i]+dt*(3/2*Ru-1/2*Ru0)/rho;
				
				
				// VERTICAL
				ve = convective_term (method, delta, xvc[i+1], x[i-1], x[i], x[i+1], x[i+2], v[j][i-1], v[j][i], v[j][i+1], v[j][i+2]);
				vw = convective_term (method, delta, xvc[i], x[i+1], x[i], x[i-1], x[i-2], v[j][i+1], v[j][i], v[j][i-1], v[j][i-2]);
				vn = convective_term (method, delta, yvc[j+1], y[j-1], y[j], y[j+1], y[j+2], v[j-1][i], v[j][i], v[j+1][i], v[j+2][i]);
				vs = convective_term (method, delta, yvc[j], y[j+1], y[j], y[j-1], y[j-2], v[j+1][i], v[j][i], v[j-1][i], v[j-2][i]);
				
				// R (vertical)
				Rv = mu*(v[j][i+1]-v[j][i])*Sv[j]/fabs(x[i+1]-x[i])+mu*(v[j+1][i]-v[j][i])*Sh[j]/fabs(y[j+1]-y[j])-mu*(v[j][i]-v[j][i-1])*Sv[j]/fabs(x[i]-x[i-1])-mu*(v[j][i]-v[j-1][i])*Sh[i]/fabs(y[j]-y[j-1])-(mflowe*ve+mflown*vn-mfloww*vw-mflows*vs);
				Rv0 = 0;
				
				// Intermediate velocity (vertical)
				vp[j][i] = v[j][i]+dt*(3/2*Rv-1/2*Rv0)/rho;
			}
		}
		
		
		// STEP 2 !!!
		
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<M; j++)
			{
				bp[j][i] = -(rho*up[j][i+1]*Sv[j]+rho*vp[j+1][i]*Sh[i]-rho*up[j][i-1]*Sv[j]-rho*vp[j-1][i]*Sh[i])/dt;
			}
		}
		Gauss_Seidel (ap, aw, ae, as, an, bp, x, fr, delta, N, M, p);
		
		
		resta = 0;
	}
	
	
//	for(int i = 0; i<N; i++)
//	{
//		cout<<x[i]<<endl;
//	}
	output_matrix(N, M, p);
	
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
			if(i==0)
			{
				x[i-1] = 0;
			}
			else if(i==N-1)
			{
				x[i+1] = x[i]+fabs(x[i]-x[i-1])/2;
			}
			if(j==0)
			{
				y[j-1] = 0;
			}
			else if(j==M-1)
			{
				y[j+1] = y[j]+fabs(y[j]-y[j-1])/2;
			}
			ae[j][i] = Sv[j]/fabs(x[i+1]-x[i]);
			aw[j][i] = Sv[j]/fabs(x[i]-x[i-1]);
			an[j][i] = Sh[i]/fabs(y[j+1]-y[j]);
			as[j][i] = Sh[i]/fabs(y[j]-y[j-1]);
			ap[j][i] = ae[j][i]+aw[j][i]+an[j][i]+as[j][i];
		}
	}
}


double convective_term (string method, float delta, float xf, float x1, float x2, float x3, float x4, double u1, double u2, double u3, double u4)
{
	// 1 refers to node W, 2 to node P, 3 to node E, and 4 to node EE
	double u, ud, uu, uuu;
	float xd, xu, xuu;
	float resta = 1;
	double uant = 0.1;
	
	if(method=="CDS")
	{
		u = 0.5*(u1+u2);
	}
	else if(method=="UDS")
	{
		int fe;
		while(resta>delta)
		{
			if(uant>0)
			{
				fe = 0;
			}
			else
			{
				fe = 1;
			}
			u = u1+fe*(u2-u1);
			resta = fabs(u-uant);
			uant = u;
		}
	}
	else if(method=="QUICK")
	{
		float g1, g2;
		while(resta>delta)
		{
			if(uant>0)
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
			resta = fabs(u-uant);
			uant = u;
		}
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
