#include<iostream>
#include<complex>
#include<math.h>
#include<vector>
#include<fstream>

using namespace std;


complex<double> diffusive(int k, int N, double Re, vector<complex<double> > u, bool LES, float CK);
complex<double> convective(int k, int N, vector<complex<double> > u);



int main()
{
	const int N = 20;
	const double Re = 40; // Reynolds number
	bool LES = 1; // 1 is LES, 0 is DNS
	double F = 0; // Source term (in Fourier space)
	
	double delta = 1e-6; // Precision of the simulation
	float CK = 0.05; // Kolgomorov constant
	float C1 = 0.02;
	double dt = C1*Re/pow(N,2); // Increment of time
	
	vector<complex<double> > u(N);
	vector<complex<double> > u0(N);
	
	for(double k = 0; k<N; k++)
	{
		u0[k] = 1/(k+1); // u at n
		u[k] = u0[k]; // u at n+1
	}
	
	complex<double> resta;
	double MAX = 1;
	
	
	double t = 0;
	
	while(MAX>delta)
	{
		t = t+dt;
		
		for(int k = 1; k<N; k++)
		{
			u[k] = u0[k]+(diffusive(k, N, Re, u0, LES, CK)-convective(k, N, u0)+F)*dt;
		}
		
		// Comprovation
		MAX = 0;
		for(int k = 1; k<N; k++)
		{
			resta = (u[k]-u0[k])/dt;
			if(abs(resta)>MAX)
			{
				MAX = abs(resta);
			}
		}
		
		for(int k = 1; k<N; k++)
		{
			u0[k] = u[k];
		}
	}
	cout<<"Steady state reached at t="<<t;
	
	vector<double> E(N);
	for(int k = 0; k<N; k++)
	{
		E[k] = abs(u[k]*conj(u[k]));
	}
	
	ofstream results;
    results.open("Results.dat");
    for(int k = 0; k<N; k++)
    {
    	results<<k+1<<"	"<<E[k]<<endl;
	}
    results.close();
	
	return 0;
}



// Calculation of the diffusive term
complex<double> diffusive(int k, int N, double Re, vector<complex<double> > u, bool LES, float CK)
{
	if(!LES)
	{
		return -(double(k)+1)*(double(k)+1)*u[k]/Re;
	}
	else
	{
		int m = 2; // Slope of the energy spectrum
		
		double viscosity;
		double eddy; // Eddy-viscosity
		double vinf;
		double vnon;
		double EkN = abs(u[N-1]*conj(u[N-1])); //Energy at the cutoff frequency
		
		vinf = 0.31*(5-m)*sqrt(3-m)*pow(CK,-3/2)/(m+1);
		vnon = 1+34.5*exp(-3.03*N/k);
		eddy = vinf*sqrt(EkN/N)*vnon;
		viscosity = 1/Re+eddy;
		return -(double(k)+1)*(double(k)+1)*u[k]*viscosity;
	}
}


// Calculation of the convective term
complex<double> convective(int k, int N, vector<complex<double> > u)
{
	complex<double> conv (0,0);
	complex<double> i(0,1);
	
	for(int p = -N; p<=N; p++)
	{
		int q = k+1-p;
		if(q>=-N && q<=N)
		{
			int qu = q;
			int pu = p;
			
			if(qu==0 || pu==0){}
			else if(qu<0)
			{
				qu = -q;
				conv = conv+u[pu-1]*i*double(q)*conj(u[qu-1]);
			}
			else if(pu<0)
			{
				pu = -p;
				conv = conv+conj(u[pu-1])*i*double(q)*u[qu-1];
			}
			else
			{
				conv = conv+u[pu-1]*i*double(q)*u[qu-1];
			}
		}
	}
	return conv;
}
