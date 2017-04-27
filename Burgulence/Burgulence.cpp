#include<iostream>
#include<complex>
#include<vector>
#include<fstream>

using namespace std;


complex<double> diffusive(int k, double Re, complex<double> u);
complex<double> convective(int k, int N, vector<complex<double> > u);
complex<double> function(int k, int N, double Re, double F, vector<complex<double> > u);



int main()
{
	const int N = 20;
	const double Re = 40;
	double F = 0;
	
	double delta = 0.0001;
	double dt = 0.005*Re/pow(N,2);
	
	vector<complex<double> > u(N);
	vector<complex<double> > u0(N);
	vector<complex<double> > uant(N);
	
	for(double k = 0; k<N; k++)
	{
		u0[k] = 1/(k+1); // u at n
		u[k] = u0[k]; // u at n+1
		uant[k] = u[k]; // calculated u at n+1
	}
	
	complex<double> resta;
	double MAX = 1;
	double MAX2;
	
	
	double t = 0;
	
	while(MAX>delta)
	{
		t = t+dt;
		cout<<endl<<t<<":"<<endl;
		
		for(int k = 1; k<N; k++)
		{
			u[k] = u0[k]+(diffusive(k, Re, u0[k])-convective(k, N, u0)+F)*dt;
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
//		cout<<MAX;
		for(int k = 1; k<N; k++)
		{
			u0[k] = u[k];
		}
	}
	
	vector<double> E(N);
	for(int k = 0; k<N; k++)
	{
		E[k] = abs(u[k]*conj(u[k]));
//		cout<<E[k]<<endl;
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




complex<double> diffusive(int k, double Re, complex<double> u)
{
	return -pow(k+1,2)*u/Re;
}


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


complex<double> function(int k, int N, double Re, double F, vector<complex<double> > u)
{
	return diffusive(k, Re, u[k])-convective(k, N, u)+F;
}
