#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

// Dimensions
const int M1 = 40;
const int M2 = 30;
const int M3 = 10;
const int N1 = 50;
const int N2 = 60;

// Definition of types
typedef double matrix[M1+M2+M3][N1+N2];


// FUNCTIONS
void horizontal_coordinates (double dx1, double dx2, double xvc[], double x[]);
void vertical_coordinates (double dy1, double dy2, double dy3, double yvc[], double y[]);
void volume (double *xvc, double *yvc, int N, int M, matrix& V);
void surface (double *yvc, int M, double Sx[]);
void properties (double *x, double *y, const float p[3][2], const float rhod[4], const float cpd[4], const float lamd[4], matrix& rho, matrix& cp, matrix& lambda);
void harmonic_mean (matrix lambda, double* x, double* y, double* xvc, double* yvc, int N, int M, matrix& lambdaw, matrix& lambdae, matrix& lambdas, matrix& lambdan);
void search_index (float point, double *x, int Number, int &ipoint, int& ip);
void constant_coefficients (double *x, double *y, double *xvc, double *yvc, double *Sx, double *Sy, matrix V, float dt, float beta, float alpha, matrix rho, matrix cp, matrix lambda, matrix lambdaw, matrix lambdae, matrix lambdas, matrix lambdan, matrix& ap, matrix& aw, matrix& ae, matrix& as, matrix& an);
void bp_coefficients (double *x, double *y, double *xvc, double *yvc, double *Sx, double *Sy, double Sytotal, matrix V, float dt, float beta, float alpha, float Qtop, float qv, float Tbottom, float Tgleft, float Tright, float Trightant, matrix Tant, matrix rho, matrix cp, matrix lambda, matrix lambdaw, matrix lambdae, matrix lambdas, matrix lambdan, matrix& bp);
void Gauss_Seidel (matrix ap, matrix aw, matrix ae, matrix as, matrix an, matrix bp, float fr, float delta, int N, int M, matrix& T);
double double_interpolation (float x, float y, double T11, double T12, double T21, double T22, double x1, double x2, double y1, double y2);
void print_matrix (matrix T, int N, int M);
void output_file (double* Tpoint1, double* Tpoint2, int Time, float dt);


int main(){
	
	// DATA
	// Coordinates
	const float p[3][2] = {
	{0.50,0.40},
	{0.50,0.70},
	{1.10,0.80}
	}; // [m]
	
	// Physical properties
	const float rhod[4] = {1500.00,1600.00,1900.00,2500.00}; // [kg/m^3]
	const float cpd[4] = {750.00,770.00,810.00,930.00}; // [J/(kgK)]
	const float lamd[4] = {170.00,140.00,200.00,140.00}; // [W/(mK)]
	
	// Boundary conditions
	const float Tbottom = 23.00; // [C]
	const float Qtop = 60.00; // [W/m]
	const float Tgleft = 33.00; // [C]
	const float alpha = 9.00; // [W/(m^2K)]
	const float Tright0 = 8.00; // Initial temperature on the right [C]
	const float variationright = 0.005; // Variation of the temperature on the right [s/C]
	const float T0 = 8.00; // Initial temperature [C]
	const float qv = 0; // Internal heat [W/m^3]
	
	// Results (coordinates)
	const float point[2][2] = {
	{0.65,0.56},
	{0.74,0.72}
	}; // Points to be studied [m]
	
	// Mathematical properties
	const int Time = 5001; // Time discretization
	const float beta = 0.5;
	const float tfinal = 5000; // Time of the simulation
	const float delta = 0.001; // Precision of the simulation
	const float fr = 1.2; // Relaxation factor
	
	
	cout<<"Program started"<<endl;
	
	// PREVIOUS CALCULATIONS
	
	float L1,L2,H1,H2,H3; // Dimensions
	L1 = p[0][0];
	L2 = p[2][0]-L1;
	H1 = p[0][1];
	H2 = p[1][1]-H1;
	H3 = p[2][1]-H1-H2;
	
	double dx1, dx2, dy1, dy2, dy3, dt; // Increments of space and time
	dt = tfinal/(Time-1); // Increment of time
	dx1 = L1/N1; // Increments in the horizontal direction
	dx2 = L2/N2;
	dy1 = H1/M1; // Increments in the vertical direction
	dy2 = H2/M2;
	dy3 = H3/M3;
	
	// Coordinates
	double xvc[N1+N2+1],yvc[M1+M2+M3+1]; // Coordinates of the faces
	double x[N1+N2],y[M1+M2+M3]; // Coordinates of the nodes
	xvc[0] = 0;
	horizontal_coordinates (dx1, dx2, xvc, x);
	yvc[0] = p[2][1];
	vertical_coordinates (dy1, dy2, dy3, yvc, y);
	
	// Surfaces and volumes
	double Sx[M1+M2+M3], Sy[N1+N2], V[M1+M2+M3][N1+N2], Sytotal; // Surfaces and volumes
	Sytotal = p[2][1]; // Total surface of the north face
	volume (xvc, yvc, N1+N2, M1+M2+M3, V);
	surface (yvc, M1+M2+M3, Sx);
	surface (xvc, N1+N2, Sy);
	
	
	cout<<"Calculating properties..."<<endl;
	
	// Density, specific heat and conductivity
	matrix rho, cp, lambda; // Density, specific heat and conductivity
	properties (x, y, p, rhod, cpd, lamd, rho, cp, lambda);
	
	
	// Harmonic mean
	matrix lambdaw, lambdae, lambdas, lambdan; // Harmonic mean
	harmonic_mean (lambda, x, y, xvc, yvc, N1+N2, M1+M2+M3, lambdaw, lambdae, lambdas, lambdan);
	
	// INITIALIZATION
	matrix T, Tant; // Temperature and Temperature in the previous instant of time
	float Tright, Trightant; // Temperature on the right and Temperature on the right in the previous instant of time
	double Tpoint1[Time], Tpoint2[Time]; // Temperatures at the points that are going to be studied
	for(int i = 0; i<N1+N2; i++)
	{
		for(int j = 0; j<M1+M2+M3; j++)
		{
			T[j][i] = T0;
			Tant[j][i] = T0;
			Tpoint1[0] = T0;
			Tpoint2[0] = T0;
		}
	}
	Tright = Tright0;
	
	// Searching for the points (0.65, 0.56) and (0.74, 0.72)
	int ipoint1, jpoint1, ip1, jp1, ipoint2, jpoint2, ip2, jp2;
    search_index (point[0][0], x, N1+N2, ipoint1, ip1);
    search_index (point[1][0], x, N1+N2, ipoint2, ip2);
    search_index (point[0][1], y, M1+M2+M3, jpoint1, jp1);
    search_index (point[1][1], y, M1+M2+M3, jpoint2, jp2);
	
	
	// CALCULATION OF CONSTANT COEFFICIENTS
	matrix ap, ae, aw, as, an, bp; // Coefficients
	constant_coefficients (x, y, xvc, yvc, Sx, Sy, V, dt, beta, alpha, rho, cp, lambda, lambdaw, lambdae, lambdas, lambdan, ap, aw, ae, as, an);
	
	
	cout<<"Solving..."<<endl;
	
	float t = 0.00; // First time increment
	double resta;
	double MAX;
	int k = 0;
	while(t<=tfinal)
	{
		k = k+1;
		t = t+dt;
		Trightant = Tright;
		Tright = Tright0+variationright*t;
		
		// CALCULATION OF NON-CONSTANT COEFFICIENTS
		bp_coefficients (x, y, xvc, yvc, Sx, Sy, Sytotal, V, dt, beta, alpha, Qtop, qv, Tbottom, Tgleft, Tright, Trightant, Tant, rho, cp, lambda, lambdaw, lambdae, lambdas, lambdan, bp);
	
		// SOLVER
		Gauss_Seidel (ap, aw, ae, as, an, bp, fr, delta, N1+N2, M1+M2+M3, T);
		
		// Assignation of the instant of time
		for(int i = 0; i<N1+N2; i++)
		{
			for(int j = 0; j<M1+M2+M3; j++)
			{
				Tant[j][i] = T[j][i];
			}
		}
		
		// Temperature at the given points
		Tpoint1[k] = double_interpolation(point[0][0], point[0][1], T[jpoint1][ipoint1], T[jp1][ipoint1], T[jpoint1][ip1], T[jp1][ip1], x[ipoint1], x[ip1], y[jpoint1], y[jp1]);
		Tpoint2[k] = double_interpolation(point[1][0], point[1][1], T[jpoint2][ipoint2], T[jp2][ipoint2], T[jpoint2][ip2], T[jp2][ip2], x[ipoint2], x[ip2], y[jpoint2], y[jp2]);
	}
	
	cout<<endl<<endl<<"Final temperature:"<<endl;
	
	// Output of the matrix temperature at the final instant of time
	print_matrix (T, N1+N2, M1+M2+M3);

    // Output file
    cout<<"Creating file..."<<endl;
    output_file (Tpoint1, Tpoint2, Time, dt);
//	resultaats (x, y, T, N1+N2, M1+M2+M3);
    
    cout<<"End of program"<<endl;
    
    ofstream results;
    results.open("Ressultats5000.dat");
    int N = N1+N2;
    int M = M1+M2+M3;
    for(int i = -1; i<N+1; i++)
    {
    	for(int j = -1; j<M+1; j++)
    	{
    		if(i==-1 && j==-1)
    		{
    			results<<0.000<<"	"<<0.800<<"	"<<(200*T[0][0]/0.005+alpha*Tgleft)/(alpha+200/0.005)<<endl;
			}
			else if(i==-1 && j==M)
			{
				results<<0.000<<"	"<<0.000<<"	"<<23.000<<endl;
			}
			else if(i==-1 && j!=-1 && j!=M)
			{
				results<<0.000<<"	"<<y[j]<<"	"<<(lambda[j][0]*T[j][0]/0.005+alpha*Tgleft)/(alpha+lambda[j][0]/0.005)<<endl;
			}
			else if(i==N && j==-1)
			{
				results<<1.100<<"	"<<0.800<<"	"<<8+0.005*tfinal<<endl;
			}
			else if(i==N && j==M)
			{
				results<<1.100<<"	"<<0.000<<"	"<<8+0.005*tfinal<<endl;
			}
			else if(i==N && j!=-1 && j!=M)
			{
				results<<1.100<<"	"<<y[j]<<"	"<<8+0.005*tfinal<<endl;
			}
			else if(j==-1 && i!=-1 && i!=N)
			{
				results<<x[i]<<"	"<<0.800<<"	"<<T[0][i]+Qtop*0.005/(1.10*lambda[0][i]*0.005)<<endl;
			}
			else if(j==M && i!=-1 && i!=N)
			{
				results<<x[i]<<"	"<<0.000<<"	"<<23.000<<endl;
			}
    		else
    		{
    			results<<x[i]<<"	"<<y[j]<<"	"<<T[j][i]<<endl;
			}
		}
		results<<endl;
	}
    results.close();
    
    return 0;
    
}



void horizontal_coordinates (double dx1, double dx2, double xvc[], double x[])
{
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
}


void vertical_coordinates (double dy1, double dy2, double dy3, double yvc[], double y[])
{
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
}


void volume (double *xvc, double *yvc, int N, int M, matrix& V)
{
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			V[j][i] = fabs(xvc[i+1]-xvc[i])*fabs(yvc[j]-yvc[j+1]); // Volume
		}
	}
}


void surface (double *yvc, int M, double Sx[])
{
	for(int j = 0; j<M; j++)
	{
		Sx[j] = fabs(yvc[j]-yvc[j+1]);
	}
}


void properties (double *x, double *y, const float p[3][2], const float rhod[4], const float cpd[4], const float lamd[4], matrix& rho, matrix& cp, matrix& lambda)
{
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
}


void harmonic_mean (matrix lambda, double* x, double* y, double* xvc, double* yvc, int N, int M, matrix& lambdaw, matrix& lambdae, matrix& lambdas, matrix& lambdan)
{
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
}


// Searching the index of the node closest to a given point (and the second closest)
void search_index (float point, double *x, int Number, int &ipoint, int& ip)
{
	for(int i = 0; i<Number-1; i++)
    {
    	if(x[i+1]-x[i]>0)
    	{
    		if(x[i]<=point && x[i+1]>point)
    		{
    			if(point-x[i]<x[i+1]-point)
				{
					ipoint = i; //ipoint is the index of the node closest to the point we want
					ip = i+1; //ip is the second node closest to it (used in interpolation)
				}
				else
				{
					ipoint = i+1;
					ip = i;
				}
			}    		
		}
		else
		{
			if(x[i]>point && x[i+1]<=point)
    		{
    			if(point-x[i+1]<x[i]-point)
    			{
    				ipoint = i;
    				ip = i+1;
				}
				else
				{
					ipoint = i+1;
					ip = i;
				}
			}
		}
	}
	
}


// Calculation of the constant coefficients
void constant_coefficients (double *x, double *y, double *xvc, double *yvc, double *Sx, double *Sy, matrix V, float dt, float beta, float alpha, matrix rho, matrix cp, matrix lambda, matrix lambdaw, matrix lambdae, matrix lambdas, matrix lambdan, matrix& ap, matrix& aw, matrix& ae, matrix& as, matrix& an)
{
	for(int i =0; i<N1+N2; i++)
	{
		for(int j = 0; j<M1+M2+M3; j++)
		{
			if(i==0 && j==0)
			{
				ae[j][i] = beta*lambdae[j][i]*Sx[j]/(x[i+1]-x[i]);
				aw[j][i] = 0;
				as[j][i] = beta*lambdas[j][i]*Sy[i]/(y[j]-y[j+1]);
				an[j][i] = 0;
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
			}
			else if(i==0 && j!=0 && j!=M1+M2+M3-1)
			{
				ae[j][i] = beta*lambdae[j][i]*Sx[j]/(x[i+1]-x[i]);
				aw[j][i] = 0;
				as[j][i] = beta*lambdas[j][i]*Sy[i]/(y[j]-y[j+1]);
				an[j][i] = beta*lambdan[j][i]*Sy[i]/(y[j-1]-y[j]);
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i]);
			}
			else if(i==0 && j==M1+M2+M3-1)
			{
				ae[j][i] = beta*lambdae[j][i]*Sx[j]/(x[i+1]-x[i]);
				aw[j][i] = 0;
				as[j][i] = 0;
				an[j][i] = beta*lambdan[j][i]*Sy[i]/(y[j-1]-y[j]);
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+beta*lambda[j][i]/(y[j]-yvc[j+1])*Sy[i];
			}
			else if(i==N1+N2-1 && j==0)
			{
				ae[j][i] = 0;
				aw[j][i] = beta*lambdaw[j][i]*Sx[j]/(x[i]-x[i-1]);
				as[j][i] = beta*lambdas[j][i]*Sy[i]/(y[j]-y[j+1]);
				an[j][i] = 0;
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j]/(xvc[i+1]-x[i]);
			}
			else if(i==N1+N2-1 && j==M1+M2+M3-1)
			{
				ae[j][i] = 0;
				aw[j][i] = beta*lambdaw[j][i]*Sx[j]/(x[i]-x[i-1]);
				as[j][i] = 0;
				an[j][i] = beta*lambdan[j][i]*Sy[i]/(y[j-1]-y[j]);
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j]/(xvc[i+1]-x[i])+beta*lambda[j][i]/(y[j]-yvc[j+1])*Sy[i];
			}
			else if(i==N1+N2-1 && j!=0 && j!=M1+M2+M3-1)
			{
				ae[j][i] = 0;
				aw[j][i] = beta*lambdaw[j][i]*Sx[j]/(x[i]-x[i-1]);
				as[j][i] = beta*lambdas[j][i]*Sy[i]/(y[j]-y[j+1]);
				an[j][i] = beta*lambdan[j][i]*Sy[i]/(y[j-1]-y[j]);
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sx[j]/(xvc[i+1]-x[i]);
			}
			else if(i!=0 && i!=N1+N2-1 && j==0)
			{
				ae[j][i] = beta*lambdae[j][i]*Sx[j]/(x[i+1]-x[i]);
				aw[j][i] = beta*lambdaw[j][i]*Sx[j]/(x[i]-x[i-1]);
				as[j][i] = beta*lambdas[j][i]*Sy[i]/(y[j]-y[j+1]);
				an[j][i] = 0;
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
			}
			else if(i!=0 && i!=N1+N2-1 && j==M1+M2+M3-1)
			{
				ae[j][i] = beta*lambdae[j][i]*Sx[j]/(x[i+1]-x[i]);
				aw[j][i] = beta*lambdaw[j][i]*Sx[j]/(x[i]-x[i-1]);
				as[j][i] = 0;
				an[j][i] = beta*lambdan[j][i]*Sy[i]/(y[j-1]-y[j]);
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt+beta*lambda[j][i]*Sy[i]/(y[j]-yvc[j+1]);
			}
			else
			{
				ae[j][i] = beta*lambdae[j][i]*Sx[j]/(x[i+1]-x[i]);
				aw[j][i] = beta*lambdaw[j][i]*Sx[j]/(x[i]-x[i-1]);
				as[j][i] = beta*lambdas[j][i]*Sy[i]/(y[j]-y[j+1]);
				an[j][i] = beta*lambdan[j][i]*Sy[i]/(y[j-1]-y[j]);
				ap[j][i] = ae[j][i]+aw[j][i]+as[j][i]+an[j][i]+rho[j][i]*cp[j][i]*V[j][i]/dt;
			}
		}
	}
}


// Calculation of non-constant coefficients
void bp_coefficients (double *x, double *y, double *xvc, double *yvc, double *Sx, double *Sy, double Sytotal, matrix V, float dt, float beta, float alpha, float Qtop, float qv, float Tbottom, float Tgleft, float Tright, float Trightant, matrix Tant, matrix rho, matrix cp, matrix lambda, matrix lambdaw, matrix lambdae, matrix lambdas, matrix lambdan, matrix& bp)
{
	for(int i =0; i<N1+N2; i++)
		{
			for(int j = 0; j<M1+M2+M3; j++)
			{
				if(i==0 && j==0)
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*((Tgleft-Tant[j][i])*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(Tant[j][i+1]-Tant[j][i])*Sx[j]/(x[i+1]-x[i])+lambdas[j][i]*(Tant[j+1][i]-Tant[j][i])*Sy[i]/(y[j]-y[j+1]))+beta*Tgleft*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+Qtop*Sy[i]/Sytotal+qv*V[j][i];
				}
				else if(i==0 && j!=0 && j!=M1+M2+M3-1)
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*((Tgleft-Tant[j][i])*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(Tant[j][i+1]-Tant[j][i])*Sx[j]/(x[i+1]-x[i])+lambdas[j][i]*(Tant[j+1][i]-Tant[j][i])*Sy[i]/(y[j]-y[j+1])+lambdan[j][i]*(Tant[j-1][i]-Tant[j][i])*Sy[i]/(y[j-1]-y[j]))+beta*Tgleft*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+qv*V[j][i];
				}
				else if(i==0 && j==M1+M2+M3-1)
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*((Tgleft-Tant[j][i])*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+lambdae[j][i]*(Tant[j][i+1]-Tant[j][i])*Sx[j]/(x[i+1]-x[i])+lambda[j][i]*(Tbottom-Tant[j][i])/(y[j]-yvc[j+1])*Sy[i]+lambdan[j][i]*(Tant[j-1][i]-Tant[j][i])*Sy[i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tbottom/(y[j]-yvc[j+1])*Sy[i]+beta*Tgleft*Sx[j]/(1/alpha+(x[i]-xvc[i])/lambda[j][i])+qv*V[j][i];
				}
				else if(i==N1+N2-1 && j==0)
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(Tant[j][i-1]-Tant[j][i])*Sx[j]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-Tant[j][i])*Sx[j]/(xvc[i+1]-x[i])+lambdas[j][i]*(Tant[j+1][i]-Tant[j][i])*Sy[i]/(y[j]-y[j+1]))+Qtop*Sy[i]/Sytotal+beta*lambda[j][i]*Tright*Sx[j]/(xvc[i+1]-x[i])+qv*V[j][i];
				}
				else if(i==N1+N2-1 && j==M1+M2+M3-1)
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(Tant[j][i-1]-Tant[j][i])*Sx[j]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-Tant[j][i])*Sx[j]/(xvc[i+1]-x[i])+lambda[j][i]*(Tbottom-Tant[j][i])/(y[j]-yvc[j+1])*Sy[i]+lambdan[j][i]*(Tant[j-1][i]-Tant[j][i])*Sy[i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tright*Sx[j]/(xvc[i+1]-x[i])+beta*lambda[j][i]*Tbottom/(y[j]-yvc[j+1])*Sy[i]+qv*V[j][i];
				}
				else if(i==N1+N2-1 && j!=0 && j!=M1+M2+M3-1)
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(Tant[j][i-1]-Tant[j][i])*Sx[j]/(x[i]-x[i-1])+lambda[j][i]*(Trightant-Tant[j][i])*Sx[j]/(xvc[i+1]-x[i])+lambdas[j][i]*(Tant[j+1][i]-Tant[j][i])*Sy[i]/(y[j]-y[j+1])+lambdan[j][i]*(Tant[j-1][i]-Tant[j][i])*Sy[i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tright*Sx[j]/(xvc[i+1]-x[i])+qv*V[j][i];
				}
				else if(i!=0 && i!=N1+N2-1 && j==0)
				{
				bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(Tant[j][i-1]-Tant[j][i])*Sx[j]/(x[i]-x[i-1])+lambdae[j][i]*(Tant[j][i+1]-Tant[j][i])*Sx[j]/(x[i+1]-x[i])+lambdas[j][i]*(Tant[j+1][i]-Tant[j][i])*Sy[i]/(y[j]-y[j+1]))+Qtop*Sy[i]/Sytotal+qv*V[j][i];
				}
				else if(i!=0 && i!=N1+N2-1 && j==M1+M2+M3-1)
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(Tant[j][i-1]-Tant[j][i])*Sx[j]/(x[i]-x[i-1])+lambdae[j][i]*(Tant[j][i+1]-Tant[j][i])*Sx[j]/(x[i+1]-x[i])+lambda[j][i]*(Tbottom-Tant[j][i])*Sy[i]/(y[j]-yvc[j+1])+lambdan[j][i]*(Tant[j-1][i]-Tant[j][i])*Sy[i]/(y[j-1]-y[j]))+beta*lambda[j][i]*Tbottom*Sy[i]/(y[j]-yvc[j+1])+qv*V[j][i];
				}
				else
				{
					bp[j][i] = rho[j][i]*cp[j][i]*Tant[j][i]*V[j][i]/dt+(1-beta)*(lambdaw[j][i]*(Tant[j][i-1]-Tant[j][i])*Sx[j]/(x[i]-x[i-1])+lambdae[j][i]*(Tant[j][i+1]-Tant[j][i])*Sx[j]/(x[i+1]-x[i])+lambdas[j][i]*(Tant[j+1][i]-Tant[j][i])*Sy[i]/(y[j]-y[j+1])+lambdan[j][i]*(Tant[j-1][i]-Tant[j][i])*Sy[i]/(y[j-1]-y[j]))+qv*V[j][i];
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


// Double interpolation
double double_interpolation (float x, float y, double T11, double T12, double T21, double T22, double x1, double x2, double y1, double y2)
{
	double result1, result2, finalresult;
	result1 = T11+(T21-T11)*(x-x1)/(x2-x1);
	result2 = T12+(T22-T12)*(x-x1)/(x2-x1);
	finalresult = result1+(result2-result1)*(y-y1)/(y2-y1);
	return finalresult;
}


// Print matrix
void print_matrix (matrix T, int N, int M)
{
	for(int j = 0; j<M; j++)
    {
        for(int i = 0; i<N; i++)
        {
            cout<<T[j][i]<<"	";  // display the current element out of the array
        }
		cout<<endl;  // go to a new line
    }
}

// Create an output file with the results
void output_file (double* Tpoint1, double* Tpoint2, int Time, float dt)
{
	ofstream puntss;
    puntss.open("Punts.dat");
    float t = 0;
    for(int k = 0; k<Time; k++)
    {
    	puntss<<t<<"	"<<Tpoint1[k]<<"	"<<Tpoint2[k]<<"\n";
    	t = t+dt;
	}
    puntss.close();
}
