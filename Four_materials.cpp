#include<iostream>
#include<math.h>

using namespace std;

int main(){
	
	
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
	int M = 10; // Vertical discretization
	int N = 10; // Horizontal discretization
	float beta = 0.5;
	float tfinal = 5000; // Time of the simulation
	float dt = 10; // Increment of time
	float delta = 0.001; // Precision of the simulation
	
	
	
	// PREVIOUS CALCULATIONS
	float deltax, deltay;
	deltax = p[2][0]/N;
	deltay = p[2][1]/M;
	
	// Coordinates
	float xvc[N+1],yvc[M+1]; // Coordinates of the faces
	float x[N],y[M]; // Coordinates of the nodes
	
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
	float Sx[M][N], Sy[M][N], V[M][N];
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			V[j][i] = (xvc[i+1]-xvc[i])*(yvc[j]-yvc[j+1]); // Volume
			Sx[j][i] = yvc[j]-yvc[j+1]; // Surfaces east and west
			Sy[j][i] = xvc[i+1]-xvc[i]; // Surfaces north and south
		}
	}
	
	// Density, specific heat and conductivity
	float rho[M][N],cp[M][N],lam[M][N];
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			if(x[i]<=p[0][0] && y[j]<=p[0][1])
			{
				rho[j][i] = rhod[0];
				cp[j][i] = cpd[0];
				lam[j][i] = lamd[0];
			}
			else if(x[i]<=p[0][0] && y[j]>p[0][1])
			{
				rho[j][i] = rhod[2];
				cp[j][i] = cpd[2];
				lam[j][i] = lamd[2];
			}
			else if(x[i]>p[0][0] && y[j]<=p[1][1])
			{
				rho[j][i] = rhod[1];
				cp[j][i] = cpd[1];
				lam[j][i] = lamd[1];
			}
			else
			{
				rho[j][i] = rhod[3];
				cp[j][i] = cpd[3];
				lam[j][i] = lamd[3];
			}
		}
	}
	float lambda[M][N]; // falta calcular les mitjanes harmòniques
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			lambda[j][i] = lam[j][i];
		}
	}
	
	
	// INITIALIZATION
	float T[M][N];
	float Tant[M][N];
	float Trightant = Tright;
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			T[j][i] = T0;
			Tant[j][i] = T[j][i];
		}
	}
	
	
	// SOLVER
	float ap[M][N],ae[M][N],aw[M][N],as[M][N],an[M][N],bp[M][N];
	float t = dt; // First time increment
	float resta = 1;
	
	while(t<=tfinal)
	{
		Tright = 8.00+0.005*t;
		while(resta>delta)
		{
			for(int i =0; i<N; i++)
			{
				for(int j = 0; j<M; j++)
				{
					if(i==0 && j==0)
					{
						ae[j][i] = beta*lambda[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = beta*lambda[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
						bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(-Tant[j][i]*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambda[j][i+1]*(Tant[j][i+1]-Tant[j][i])/(x[i+1]-x[i])+lambda[j+1][i]*(Tant[j+1][i]-Tant[j][i])/(y[j]-y[j+1]))+Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+Qtop*Sy[j][i]/p[2][1];
					}
					else if(i==0 && j!=0 && j!=M)
					{
						ae[j][i] = beta*lambda[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = beta*lambda[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambda[j-1][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
						bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(-Tant[j][i]*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambda[j][i+1]*(Tant[j][i+1]-Tant[j][i])/(x[i+1]-x[i])+lambda[j+1][i]*(Tant[j+1][i]-Tant[j][i])/(y[j]-y[j+1])+lambda[j-1][i]*(Tant[j-1][i]-Tant[j][i])/(y[j-1]-y[j]))+Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
					}
					else if(i==0 && j==M)
					{
						ae[j][i] = beta*lambda[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = 0;
						as[j][i] = 0;
						an[j][i] = beta*lambda[j-1][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+beta*lambda[j][i]/(y[j]-yvc[j+1])*Sy[j][i];
						bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(-Tant[j][i]*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambda[j][i+1]*(Tant[j][i+1]-Tant[j][i])/(x[i+1]-x[i])-lambda[j][i]*Tant[j][i]/(y[j]-yvc[j+1])*Sy[j][i]+lambda[j-1][i]*(Tant[j-1][i]-Tant[j][i])/(y[j-1]-y[j]))+lambda[j][i]*Tbottom/(y[j]-yvc[j+1])*Sy[j][i]+Tgleft*Sx[j][i]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
					}
					else if(i!=0 && i!=N && j==0)
					{
						ae[j][i] = beta*lambda[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambda[j][i-1]/(x[i]-x[i-1]);
						as[j][i] = beta*lambda[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
						bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*((lambda[j][i-1]*(Tant[j][i-1]-Tant[j][i])/(x[i]-x[i-1])+(x[i]-xvc[i])/lambda[j][i])+lambda[j][i+1]*(Tant[j][i+1]-Tant[j][i])/(x[i+1]-x[i])+lambda[j+1][i]*(Tant[j][i]-Tant[j+1][i])/(y[j]-y[j+1]))+Qtop*Sy[j][i]/p[2][1];
					}
					else if(i==N && j==0)
					{
						ae[j][i] = 0;
						aw[j][i] = beta*lambda[j][i-1]/(x[i]-x[i-1]);
						as[j][i] = beta*lambda[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = 0;
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j][i+1]/(xvc[i+1]-x[i]);
						bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*((lambda[j][i-1]*(Tant[j][i-1]-Tant[j][i])/(x[i]-x[i-1])+(x[i]-xvc[i])/lambda[j][i])+lambda[j][i]*(Trightant-Tant[j][i])*Sx[j][i+1]/(xvc[i+1]-x[i])+lambda[j+1][i]*(Tant[j][i]-Tant[j+1][i])/(y[j]-y[j+1]))+Qtop*Sy[j][i]/p[2][1]+lambda[j][i]*Tright*Sx[j][i+1]/(xvc[i+1]-x[i]);
					}
					// falta fer uns quants casos...
					else
					{
						ae[j][i] = beta*lambda[j][i+1]/(x[i+1]-x[i]);
						aw[j][i] = beta*lambda[j][i-1]/(x[i]-x[i-1]);
						as[j][i] = beta*lambda[j+1][i]/(y[j]-y[j+1]);
						an[j][i] = beta*lambda[j-1][i]/(y[j-1]-y[j]);
						ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
						bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(lambda[j][i-1]*(Tant[j][i-1]-Tant[j][i])/(x[i]-x[i-1])+lambda[j][i+1]*(Tant[j][i+1]-Tant[j][i])/(x[i+1]-x[i])+lambda[j+1][i]*(Tant[j+1][i]-Tant[j][i])/(y[j]-y[j+1])+lambda[j-1][i]*(Tant[j-1][i]-Tant[j][i])/(y[j-1]-y[j]));
					}
				}
			}
			
			resta = delta;
		}
		
		t = t+dt;
	}
	cout<<"\n\n";
	
	
	
	// SCREEN ! ! ! ! ! ! ! ! :D
	// showing the matrix on the screen
    for(int j=0;j<M;j++)  // loop 3 times for three lines
    {
        for(int i=0;i<N;i++)  // loop for the three elements on the line
        {
            cout<<lambda[j][i]<<" ";  // display the current element out of the array
        }
    cout<<endl;  // when the inner loop is done, go to a new line
    }
    
    cout<<"\nFighting!~ ^w^";
}
