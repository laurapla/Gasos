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
	int N = 20;
	double Re = 40;
	double F = 0;
	
	double delta = 0.00001;
	double dt = 0.1*Re/pow(N,2);
	
	vector<complex<double> > u(N);
	vector<complex<double> > u0(N);
	vector<complex<double> > uant(N);
	
	for(double k = 1; k<N; k++)
	{
		u0[k] = 1/k;
		u[k] = u0[k];
		uant[k] = u[k];
	}
	
	complex<double> resta;
	double MAX = 1;
	double MAX2;
	
	complex<double> diff;
	complex<double> conv;
	
	
	double t = 0;
	
	while(MAX>delta)
	{
		t = t+dt;
		cout<<endl<<t<<":"<<endl;
		
		for(int k = 1; k<N; k++)
		{
			u[k] = u0[k]+(1.5*function(k, N, Re, F, u)-0.5*function(k, N, Re, F, u0));
		}
		
		// Comprovation
		MAX = 0;
		for(int k = 0; k<N; k++)
		{
			resta = u[k]-u0[k];
			if(abs(resta)>MAX)
			{
				MAX = abs(resta);
			}
		}
		cout<<MAX;
		for(int k = 0; k<N; k++)
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
    for(int k = 1; k<N; k++)
    {
    	results<<k<<"	"<<E[k]<<endl;
	}
    results.close();
	
	return 0;
}




complex<double> diffusive(int k, double Re, complex<double> u)
{
	return -pow(k,2)*u/Re;
}


complex<double> convective(int k, int N, vector<complex<double> > u)
{
	complex<double> conv (0,0);
	complex<double> i(0,1);
	
	for(int p = -N+1; p<N; p++)
	{
		int q = k-p;
		if(q>-N && q<N)
		{
			int qu = q;
			int pu = p;
			if(qu<0)
			{
				qu = -q;
			}
			if(pu<0)
			{
				pu = -p;
			}
			conv = conv+u[pu]*i*double(q)*u[qu];
		}
	}
	return conv;
}


complex<double> function(int k, int N, double Re, double F, vector<complex<double> > u)
{
	return diffusive(k, Re, u[k])-convective(k, N, u)+F;
}
