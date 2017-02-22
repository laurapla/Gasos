#include<iostream>

using namespace std;


// DATA
// Coordinates
float p[3][2] = {
{0.50,0.40},
{0.50,0.70},
{1.10,0.80}
}; //m

// Physical properties
float rhod[4] = {1500.00,1600.00,1900.00,2500.00}; //kg/m^3
float cpd[4] = {750.00,770.00,810.00,930.00}; //J/(kgK)
float lamd[4] = {170.00,140.00,200.00,140.00}; //W/(mK)

// Boundary conditions
float Tbottom = 23.00; //ºC
float Qtop = 60.00; //W/m
float Tgleft = 33.00; //ºC
float alpha = 9.00; //W/(m^2K)
float Tright = 8.00; //ºC
float T0 = 8.00; //ºC Initial temperature

// Mathematical properties
const int M = 10; // Vertical discretization
const int N = 10; // Horizontal discretization
const int L = 10001; // Time discretization
float beta = 0.5;
float tfinal = 10000; // Time of the simulation
float delta = 0.001; // Precision of the simulation
float fr = 1.2; // Relaxation factor


float deltax, deltay, dt; // Increments of space and time
float xvc[N+1],yvc[M+1]; // Coordinates of the faces
float x[N],y[M]; // Coordinates of the nodes
double Sx[M][N+1], Sy[M+1][N], V[M][N]; // Surfaces and volumes
float rho[M][N],cp[M][N],lambda[M][N]; // Density, specific heat and conductivity
double lambdaw[M][N], lambdae[M][N], lambdas[M][N], lambdan[M][N]; // Harmonic mean
double T[L][M][N]; // Temperature
double Tant[M][N]; // Temperature on the previous instant of time
double Tcalc[M][N]; // Calculated temperature
float Trightant; // Temperature on the right in the previous instant of time
double ap[M][N],ae[M][N],aw[M][N],as[M][N],an[M][N],bp[M][N]; // Coefficients
float t = 0.00; // First time increment
double resta;
double MAX;
int k = 0;
float point[2][2] = {
{0.65,0.56},
{0.74,0.72}
};
int ipoint1, jpoint1, ipoint2, jpoint2;



int main(){
	
	
	// PREVIOUS CALCULATIONS
	
	dt = tfinal/(L-1); // Increment of time
	deltax = p[2][0]/N;
	deltay = p[2][1]/M;
	
	// Coordinates	
	xvc[0] = 0.00;
	for (int i = 1; i<N+1; i++)
	{
		xvc[i] = xvc[i-1]+deltax;
		x[i-1] = (xvc[i-1]+xvc[i])/2;
	}
	
	yvc[0] = p[2][1];
	for (int j = 1; j<M+1; j++)
	{
		yvc[j] = yvc[j-1]-deltay;
		y[j-1] = (yvc[j-1]+yvc[j])/2;
	}
	
	// Surfaces and volumes
	float Sytotal = p[2][1];
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			V[j][i] = (xvc[i+1]-xvc[i])*(yvc[j]-yvc[j+1]); // Volume
		}
	}
	for(int i = 0; i<=N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			Sx[j][i] = yvc[j]-yvc[j+1]; // Surfaces east and west
		}
	}
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<=M; j++)
		{
			Sy[j][i] = xvc[i+1]-xvc[i]; // Surfaces north and south
		}
	}
	
	// Density, specific heat and conductivity
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(x[i]<=p[0][0] && y[j]<=p[0][1])
			{
				rho[j][i] = rhod[0];
				cp[j][i] = cpd[0];
				lambda[j][i] = lamd[0];
			}
			else if(x[i]<=p[0][0] && y[j]>p[0][1])
			{
				rho[j][i] = rhod[2];
				cp[j][i] = cpd[2];
				lambda[j][i] = lamd[2];
			}
			else if(x[i]>p[0][0] && y[j]<=p[1][1])
			{
				rho[j][i] = rhod[1];
				cp[j][i] = cpd[1];
				lambda[j][i] = lamd[1];
			}
			else
			{
				rho[j][i] = rhod[3];
				cp[j][i] = cpd[3];
				lambda[j][i] = lamd[3];
			}
		}
	}
	
	// Harmonic mean
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(i==0)
			{
				lambdaw[j][i] = lambda[j][i];
				lambdan[j][i] = (y[j-1]-y[j])/((y[j-1]-yvc[j])/lambda[j-1][i]+(yvc[j]-y[j])/lambda[j][i]);
				lambdae[j][i] = (x[i+1]-x[i])/((x[i+1]-xvc[i+1])/lambda[j][i+1]+(xvc[i+1]-x[i])/lambda[j][i]);
				lambdas[j][i] = (y[j]-y[j+1])/((yvc[j+1]-y[j+1])/lambda[j+1][i]+(y[j]-yvc[j+1])/lambda[j][i]);
			}
			else if(i==N-1)
			{
				lambdaw[j][i] = (x[i]-x[i-1])/((x[i]-xvc[i])/lambda[j][i-1]+(x[i]-xvc[i])/lambda[j][i]);
				lambdan[j][i] = (y[j-1]-y[j])/((y[j-1]-yvc[j])/lambda[j-1][i]+(yvc[j]-y[j])/lambda[j][i]);
				lambdae[j][i] = lambda[j][i];
				lambdas[j][i] = (y[j]-y[j+1])/((yvc[j+1]-y[j+1])/lambda[j+1][i]+(y[j]-yvc[j+1])/lambda[j][i]);
			}
			else if(j==0)
			{
				lambdaw[j][i] = (x[i]-x[i-1])/((x[i]-xvc[i])/lambda[j][i-1]+(x[i]-xvc[i])/lambda[j][i]);
				lambdan[j][i] = lambda[j][i];
				lambdae[j][i] = (x[i+1]-x[i])/((x[i+1]-xvc[i+1])/lambda[j][i+1]+(xvc[i+1]-x[i])/lambda[j][i]);
				lambdas[j][i] = (y[j]-y[j+1])/((yvc[j+1]-y[j+1])/lambda[j+1][i]+(y[j]-yvc[j+1])/lambda[j][i]);
			}
			else if(j==M-1)
			{
				lambdaw[j][i] = (x[i]-x[i-1])/((x[i]-xvc[i])/lambda[j][i-1]+(x[i]-xvc[i])/lambda[j][i]);
				lambdan[j][i] = (y[j-1]-y[j])/((y[j-1]-yvc[j])/lambda[j-1][i]+(yvc[j]-y[j])/lambda[j][i]);
				lambdae[j][i] = (x[i+1]-x[i])/((x[i+1]-xvc[i+1])/lambda[j][i+1]+(xvc[i+1]-x[i])/lambda[j][i]);
				lambdas[j][i] = lambda[j][i];
			}
			else
			{
				lambdaw[j][i] = (x[i]-x[i-1])/((x[i]-xvc[i])/lambda[j][i-1]+(x[i]-xvc[i])/lambda[j][i]);
				lambdan[j][i] = (y[j-1]-y[j])/((y[j-1]-yvc[j])/lambda[j-1][i]+(yvc[j]-y[j])/lambda[j][i]);
				lambdae[j][i] = (x[i+1]-x[i])/((x[i+1]-xvc[i+1])/lambda[j][i+1]+(xvc[i+1]-x[i])/lambda[j][i]);
				lambdas[j][i] = (y[j]-y[j+1])/((yvc[j+1]-y[j+1])/lambda[j+1][i]+(y[j]-yvc[j+1])/lambda[j][i]);
			}
		}
	}
	
	
	// INITIALIZATION
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			T[0][j][i] = T0;
			Tcalc[j][i] = T[0][j][i];
		}
	}
	
	
	while(t<=tfinal)
	{
		k = k+1;
		t = t+dt;
		Trightant = Tright;
		Tright = 8.00+0.005*t;
		
		// CALCULATION OF COEFFICIENTS
			for(int i =0; i<N; i++)
			{
				for(int j = 0; j<M; j++)
				{
					if(i==0 && j==0)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*((Tgleft-T[k-1][j][i])*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1]))+beta*Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+Qtop*Sy[j][i]/Sytotal;
					}
					else if(i==0 && j!=0 && j!=M-1)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*((Tgleft-T[k-1][j][i])*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
					}
					else if(i==0 && j==M-1)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = 0;
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+beta*lambda[j][i]/(y[j]-yvc[j+1])*Sy[j+1][i];
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*((Tgleft-T[k-1][j][i])*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambda[j][i]*(Tbottom-T[k-1][j][i])/(y[j]-yvc[j+1])*Sy[j+1][i]+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tbottom/(y[j]-yvc[j+1])*Sy[j+1][i]+beta*Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
					}
					else if(i==N-1 && j==0)
					{
						ae[j][i] = 0;
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j][i+1]/(xvc[i+1]-x[i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-T[k-1][j][i])*Sx[j][i+1]/(xvc[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1]))+Qtop*Sy[j][i]/Sytotal+beta*lambda[j][i]*Tright*Sx[j][i+1]/(xvc[i+1]-x[i]);
					}
					else if(i==N-1 && j==M-1)
					{
						ae[j][i] = 0;
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = 0;
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j][i+1]/(xvc[i+1]-x[i])+beta*lambda[j][i]/(y[j]-yvc[j+1])*Sy[j+1][i];
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-T[k-1][j][i])*Sx[j][i+1]/(xvc[i+1]-x[i])+lambda[j][i]*(Tbottom-T[k-1][j][i])/(y[j]-yvc[j+1])*Sy[j+1][i]+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tright*Sx[j][i+1]/(xvc[i+1]-x[i])+beta*lambda[j][i]*Tbottom/(y[j]-yvc[j+1])*Sy[j+1][i];
					}
					else if(i==N-1 && j!=0 && j!=M-1)
					{
						ae[j][i] = 0;
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j][i+1]/(xvc[i+1]-x[i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-T[k-1][j][i])*Sx[j][i+1]/(xvc[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tright*Sx[j][i+1]/(xvc[i+1]-x[i]);
					}
					else if(i!=0 && i!=N-1 && j==0)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1]))+Qtop*Sy[j][i]/Sytotal;
					}
					else if(i!=0 && i!=N-1 && j==M-1)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = 0;
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sy[j+1][i]/(y[j]-yvc[j+1]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambda[j][i]*(Tbottom-T[k-1][j][i])*Sy[j+1][i]/(y[j]-yvc[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tbottom*Sy[j+1][i]/(y[j]-yvc[j+1]);
					}
					else
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]));
					}
				}
			}
		
		MAX = 1;
		
		while(MAX>delta)
		{			
			
			// SOLVER: Gauss-Seidel
			for(int i = 0; i<N; i++)
			{
				for(int j = 0; j<M; j++)
				{
					if(i==0 && j==0)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==0 && j==M-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==0 && j!=0 && j!=M-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==N-1 && j==0)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==N-1 && j==M-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==N-1 && j!=0 && j!=M-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i!=0 && i!=N-1 && j==0)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i!=0 && i!=N-1 && j==M-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
				}
			}
			
			// Comprovation
			MAX = 0;
			for(int i = 0; i<N; i++)
			{
				for(int j = 0; j<M; j++)
				{
					if(Tcalc[j][i]>T[k][j][i])
					{
						resta = Tcalc[j][i]-T[k][j][i];
					}
					else
					{
						resta = T[k][j][i]-Tcalc[j][i];
					}
					
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
					Tcalc[j][i] = T[k][j][i];
				}
			}	
		}
	}
	
	
	
	// SCREEN ! ! ! ! ! ! ! ! :D
	// showing the matrix on the screen
    for(int j=0;j<M;j++)  // loop 3 times for three lines
    {
        for(int i=0;i<N;i++)  // loop for the three elements on the line
        {
            cout<<T[L][j][i]<<"	";  // display the current element out of the array
        }
		cout<<endl;  // when the inner loop is done, go to a new line
    }
    
    
    // Points
    for(int i = 0; i<N-1; i++)
    {
    	if(x[i]<=point[0][0] && x[i+1]>point[0][0])
    	{
    		if(point[0][0]-x[i]<x[i+1]-point[0][0])
    		{
    			ipoint1 = i;
			}
			else
			{
				ipoint1 = i+1;
			}
		}
		if(x[i]<=point[1][0] && x[i+1]>point[1][0])
    	{
    		if(point[1][0]-x[i]<x[i+1]-point[1][0])
    		{
    			ipoint2 = i;
			}
			else
			{
				ipoint2 = i+1;
			}
		}
	}
	for(int j = 0; j<M-1; j++)
    {
    	if(y[j]>point[0][1] && y[j+1]<=point[0][1])
    	{
    		if(point[0][1]-y[j+1]<y[j]-point[0][1])
    		{
    			jpoint1 = j;
			}
			else
			{
				jpoint1 = j+1;
			}
		}
		if(y[j]>point[1][1] && y[j+1]<=point[1][1])
    	{
    		if(point[1][1]-y[j+1]<y[j]-point[1][1])
    		{
    			jpoint2 = j;
			}
			else
			{
				jpoint2 = j+1;
			}
		}
	}
    
    // SCREEN ! ! ! ! ! ! ! ! :D
	// showing the matrix on the screen
    for(int j=0;j<1;j++)  // loop 3 times for three lines
    {
        for(int i=0;i<N;i++)  // loop for the three elements on the line
        {
            cout<<x[i]<<"	";  // display the current element out of the array
        }
		cout<<endl;  // when the inner loop is done, go to a new line
    }
    // SCREEN ! ! ! ! ! ! ! ! :D
	// showing the matrix on the screen
    for(int j=0;j<1;j++)  // loop 3 times for three lines
    {
        for(int i=0;i<M;i++)  // loop for the three elements on the line
        {
            cout<<y[i]<<"	";  // display the current element out of the array
        }
		cout<<endl;  // when the inner loop is done, go to a new line
    }
    
    cout<<"\n"<<x[ipoint1]<<","<<y[jpoint1]<<"\n"<<x[ipoint2]<<","<<y[jpoint2];
    cout<<"\nFighting!~ ^w^";
    
}



