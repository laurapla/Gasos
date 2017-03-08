#include <iostream>
#include <vector>

double mean(double *array, size_t n)
{
	double m = 0;
	for(size_t i = 0; i<n; ++i){
		m += array[i];
	}
	return m/n;
}

int main(){
	
	// version 1: array
//	double a[] = {1, 2, 3, 4, 5};
//	std::cout<<mean(a, 5)<<std::endl;	// will print 3
	
	// version 2: doing exactly the same but with vectors
//	std::vector<double> a;
//	a.push_back(1);
//	a.push_back(2);
//	a.push_back(3);
//	a.push_back(4);
//	a.push_back(5);
//	std::cout<<mean(&a[0], 5)<<std::endl;	// will print 3
	// we see that we had to push_back the elements into the vector :(
	
	// version 3: we can define the vector as a copy of the array
	double p[] = {1, 2, 3, 4, 5};
	std::vector<double> a(p, p+5);
	std::cout<<mean(&a[0], 5)<<std::endl;	// will print 3
	
	return 0;
}
