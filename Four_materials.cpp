#include<iostream>
#include<fstream>

using namespace std;


// DATA
// Coordinates
const float p[3][2] = {
{0.50,0.40},
{0.50,0.70},
{1.10,0.80}
}; // [m]
const float B = 1; // Depth [m]

// Physical properties
const float rhod[4] = {1500.00,1600.00,1900.00,2500.00}; // [kg/m^3]
const float cpd[4] = {750.00,770.00,810.00,930.00}; // [J/(kgK)]
const float lamd[4] = {170.00,140.00,200.00,140.00}; // [W/(mK)]

// Boundary conditions
const float Tbottom = 23.00; // [ºC]
const float Qtop = 60.00; // [W/m]
const float Tgleft = 33.00; // [ºC]
const float alpha = 9.00; // [W/(m^2K)]
const float Tright0 = 8.00; // Initial temperature on the right [ºC]
const float variationright = 0.005; // Variation of the temperature on the right [s/ºC]
const float T0 = 8.00; // Initial temperature [ºC]
const float qv = 0; // Internal heat [W/m^3]

// Mathematical properties
const int M1 = 20;
const int M2 = 40;
const int M3 = 20;
const int N1 = 40;
const int N2 = 60;
const int Time = 10001; // Time discretization
const float beta = 0.5;
const float tfinal = 10000; // Time of the simulation
const float delta = 0.001; // Precision of the simulation
const float fr = 1.2; // Relaxation 

// Results (coordinates)
const float point[2][2] = {
{0.65,0.56},
{0.74,0.72}
}; // Points to be studied [m]


// Declaration of variables
float L1,L2,H1,H2,H3;
double dx1, dx2, dy1, dy2, dy3, dt; // Increments of space and time
double xvc[N1+N2+1],yvc[M1+M2+M3+1]; // Coordinates of the faces
double x[N1+N2],y[M1+M2+M3]; // Coordinates of the nodes
double Sx[M1+M2+M3][N1+N2+1], Sy[M1+M2+M3+1][N1+N2], V[M1+M2+M3][N1+N2], Sytotal; // Surfaces and volumes
float rho[M1+M2+M3][N1+N2],cp[M1+M2+M3][N1+N2],lambda[M1+M2+M3][N1+N2]; // Density, specific heat and conductivity
double lambdaw[M1+M2+M3][N1+N2], lambdae[M1+M2+M3][N1+N2], lambdas[M1+M2+M3][N1+N2], lambdan[M1+M2+M3][N1+N2]; // Harmonic mean
double T[Time][M1+M2+M3][N1+N2]; // Temperature
double Tcalc[M1+M2+M3][N1+N2]; // Calculated temperature
float Tright; // Temperature on the right
float Trightant; // Temperature on the right in the previous instant of time
double ap[M1+M2+M3][N1+N2],ae[M1+M2+M3][N1+N2],aw[M1+M2+M3][N1+N2],as[M1+M2+M3][N1+N2],an[M1+M2+M3][N1+N2],bp[M1+M2+M3][N1+N2]; // Coefficients
float t = 0.00; // First time increment
double resta;
double MAX;
int k = 0;
int ipoint1, jpoint1, ipoint2, jpoint2;



int main(){
	
	cout<<"Program started"<<endl;
	
	// PREVIOUS CALCULATIONS
	L1 = p[0][0];
	L2 = p[2][0]-L1;
	H1 = p[0][1];
	H2 = p[1][1]-H1;
	H3 = p[2][1]-H1-H2;
	
	dt = tfinal/(Time-1); // Increment of time
	dx1 = L1/N1; // Increments in the horizontal direction
	dx2 = L2/N2;
	dy1 = H1/M1; // Increments in the vertical direction
	dy2 = H2/M2;
	dy3 = H3/M3;
	
	// Coordinates	
	xvc[0] = 0;
	for (int i = 1; i<N1+N2+1; i++)
	{
		if(i<=N1)
		{
			xvc[i] = xvc[i-1]+dx1;
			x[i-1] = (xvc[i-1]+xvc[i])/2;
		}
		else
		{
			xvc[i] = xvc[i-1]+dx2;
			x[i-1] = (xvc[i-1]+xvc[i])/2;
		}
	}
	
	yvc[0] = p[2][1];
	for (int j = 1; j<M1+M2+M3+1; j++)
	{
		if(j<=M3)
		{
			yvc[j] = yvc[j-1]-dy3;
			y[j-1] = (yvc[j-1]+yvc[j])/2;
		}
		else if(j>M3 && j<=M2+M3)
		{
			yvc[j] = yvc[j-1]-dy2;
			y[j-1] = (yvc[j-1]+yvc[j])/2;
		}
		else
		{
			yvc[j] = yvc[j-1]-dy1;
			y[j-1] = (yvc[j-1]+yvc[j])/2;
		}
	}
	
	// Surfaces and volumes
	Sytotal = B*p[2][1]; // Total surface of the north face
	for(int i = 0; i<N1+N2; i++)
	{
		for(int j = 0; j<M1+M2+M3; j++)
		{
			V[j][i] = B*(xvc[i+1]-xvc[i])*(yvc[j]-yvc[j+1]); // Volume
		}
	}
	for(int i = 0; i<=N1+N2; i++)
	{
		for(int j = 0; j<M1+M2+M3; j++)
		{
			Sx[j][i] = B*(yvc[j]-yvc[j+1]); // Surfaces east and west
		}
	}
	for(int i = 0; i<N1+N2; i++)
	{
		for(int j = 0; j<=M1+M2+M3; j++)
		{
			Sy[j][i] = xvc[i+1]-xvc[i]; // Surfaces north and south
		}
	}
	
	cout<<"Calculating properties..."<<endl;
	
	// Density, specific heat and conductivity
	for(int i = 0; i<N1+N2; i++)
	{
		for(int j = 0; j<M1+M2+M3; j++)
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
	for(int i = 0; i<N1+N2; i++)
	{
		for(int j = 0; j<M1+M2+M3; j++)
		{
			if(i==0)
			{
				lambdaw[j][i] = lambda[j][i];
				lambdan[j][i] = (y[j-1]-y[j])/((y[j-1]-yvc[j])/lambda[j-1][i]+(yvc[j]-y[j])/lambda[j][i]);
				lambdae[j][i] = (x[i+1]-x[i])/((x[i+1]-xvc[i+1])/lambda[j][i+1]+(xvc[i+1]-x[i])/lambda[j][i]);
				lambdas[j][i] = (y[j]-y[j+1])/((yvc[j+1]-y[j+1])/lambda[j+1][i]+(y[j]-yvc[j+1])/lambda[j][i]);
			}
			else if(i==N1+N2-1)
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
			else if(j==M1+M2+M3-1)
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
	for(int i = 0; i<N1+N2; i++)
	{
		for(int j = 0; j<M1+M2+M3; j++)
		{
			T[0][j][i] = T0;
			Tcalc[j][i] = T[0][j][i];
		}
	}
	Tright = Tright0;
	
	cout<<"Solving..."<<endl;
	
	while(t<=tfinal)
	{
		k = k+1;
		t = t+dt;
		Trightant = Tright;
		Tright = Tright0+variationright*t;
		
		// CALCULATION OF COEFFICIENTS
			for(int i =0; i<N1+N2; i++)
			{
				for(int j = 0; j<M1+M2+M3; j++)
				{
					if(i==0 && j==0)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*((Tgleft-T[k-1][j][i])*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1]))+beta*Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+Qtop*Sy[j][i]/Sytotal+qv*V[j][i];
					}
					else if(i==0 && j!=0 && j!=M1+M2+M3-1)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*((Tgleft-T[k-1][j][i])*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+qv*V[j][i];
					}
					else if(i==0 && j==M1+M2+M3-1)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = 0;
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+beta*lambda[j][i]/(y[j]-yvc[j+1])*Sy[j+1][i];
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*((Tgleft-T[k-1][j][i])*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambda[j][i]*(Tbottom-T[k-1][j][i])/(y[j]-yvc[j+1])*Sy[j+1][i]+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tbottom/(y[j]-yvc[j+1])*Sy[j+1][i]+beta*Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+qv*V[j][i];
					}
					else if(i==N1+N2-1 && j==0)
					{
						ae[j][i] = 0;
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j][i+1]/(xvc[i+1]-x[i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-T[k-1][j][i])*Sx[j][i+1]/(xvc[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1]))+Qtop*Sy[j][i]/Sytotal+beta*lambda[j][i]*Tright*Sx[j][i+1]/(xvc[i+1]-x[i])+qv*V[j][i];
					}
					else if(i==N1+N2-1 && j==M1+M2+M3-1)
					{
						ae[j][i] = 0;
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = 0;
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j][i+1]/(xvc[i+1]-x[i])+beta*lambda[j][i]/(y[j]-yvc[j+1])*Sy[j+1][i];
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-T[k-1][j][i])*Sx[j][i+1]/(xvc[i+1]-x[i])+lambda[j][i]*(Tbottom-T[k-1][j][i])/(y[j]-yvc[j+1])*Sy[j+1][i]+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tright*Sx[j][i+1]/(xvc[i+1]-x[i])+beta*lambda[j][i]*Tbottom/(y[j]-yvc[j+1])*Sy[j+1][i]+qv*V[j][i];
					}
					else if(i==N1+N2-1 && j!=0 && j!=M1+M2+M3-1)
					{
						ae[j][i] = 0;
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j][i+1]/(xvc[i+1]-x[i]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-T[k-1][j][i])*Sx[j][i+1]/(xvc[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tright*Sx[j][i+1]/(xvc[i+1]-x[i])+qv*V[j][i];
					}
					else if(i!=0 && i!=N1+N2-1 && j==0)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1]))+Qtop*Sy[j][i]/Sytotal+qv*V[j][i];
					}
					else if(i!=0 && i!=N1+N2-1 && j==M1+M2+M3-1)
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = 0;
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sy[j+1][i]/(y[j]-yvc[j+1]);
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambda[j][i]*(Tbottom-T[k-1][j][i])*Sy[j+1][i]/(y[j]-yvc[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tbottom*Sy[j+1][i]/(y[j]-yvc[j+1])+qv*V[j][i];
					}
					else
					{
						ae[j][i] = beta*lambdae[j][i]*Sx[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambdaw[j][i]*Sx[j][i]/(x[i]-x[i-1]);
						as[j][i] = beta*lambdas[j][i]*Sy[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambdan[j][i]*Sy[j][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
						bp[j][i] = rho[j][i]*cp[j][i]*T[k-1][j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(T[k-1][j][i-1]-T[k-1][j][i])*Sx[j][i]/(x[i]-x[i-1])+lambdae[j][i]*(T[k-1][j][i+1]-T[k-1][j][i])*Sx[j][i+1]/(x[i+1]-x[i])+lambdas[j][i]*(T[k-1][j+1][i]-T[k-1][j][i])*Sy[j+1][i]/(y[j]-y[j+1])+lambdan[j][i]*(T[k-1][j-1][i]-T[k-1][j][i])*Sy[j][i]/(y[j-1]-y[j]))+qv*V[j][i];
					}
				}
			}
		
		MAX = 1;
		
		while(MAX>delta)
		{			
			
			// SOLVER: Gauss-Seidel
			for(int i = 0; i<N1+N2; i++)
			{
				for(int j = 0; j<M1+M2+M3; j++)
				{
					if(i==0 && j==0)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==0 && j==M1+M2+M3-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==0 && j!=0 && j!=M1+M2+M3-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==N1+N2-1 && j==0)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==N1+N2-1 && j==M1+M2+M3-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i==N1+N2-1 && j!=0 && j!=M1+M2+M3-1)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+as[j][i]*Tcalc[j+1][i]+an[j][i]*T[k][j-1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i!=0 && i!=N1+N2-1 && j==0)
					{
						T[k][j][i] = Tcalc[j][i]+fr*((aw[j][i]*T[k][j][i-1]+ae[j][i]*Tcalc[j][i+1]+as[j][i]*Tcalc[j+1][i]+bp[j][i])/ap[j][i]-Tcalc[j][i]);
					}
					else if(i!=0 && i!=N1+N2-1 && j==M1+M2+M3-1)
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
			for(int i = 0; i<N1+N2; i++)
			{
				for(int j = 0; j<M1+M2+M3; j++)
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
			for(int i = 0; i<N1+N2; i++)
			{
				for(int j = 0; j<M1+M2+M3; j++)
				{
					Tcalc[j][i] = T[k][j][i];
				}
			}	
		}
	}
	
	cout<<endl<<endl<<"Final temperature:"<<endl;
	
	// SCREEN ! ! ! ! ! ! ! ! :D
    for(int j = 0; j<M1+M2+M3; j++)
    {
        for(int i = 0; i<N1+N2; i++)
        {
            cout<<T[Time][j][i]<<"	";  // display the current element out of the array
        }
		cout<<endl;  // go to a new line
    }
    
    
    // Points
    for(int i = 0; i<N1+N2-1; i++)
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
	for(int j = 0; j<M1+M2+M3-1; j++)
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
    
    cout<<endl<<x[ipoint1]<<","<<y[jpoint1]<<"\n"<<x[ipoint2]<<","<<y[jpoint2]<<endl;
    
    // Output file
    cout<<"Creating file..."<<endl;
    ofstream results;
    results.open("Resultats.txt");
    t = 0;
    for(int k = 0; k<Time; k++)
    {
    	results<<t<<"	"<<T[k][jpoint1][ipoint1]<<"	"<<T[k][jpoint2][ipoint2]<<"\n";
    	t = t+dt;
	}
    results.close();
    
    cout<<"End of program"<<endl;
    
    return 0;
    
}

