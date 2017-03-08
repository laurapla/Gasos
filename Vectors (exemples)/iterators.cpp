#include <iostream>
#include <vector>

int main()
{
	std::vector<double> a;
	std::vector<double>::const_iterator i;		// declares a const iterator i for a vector<double>
	a.push_back(1);
	a.push_back(2);
	a.push_back(3);
	a.push_back(4);
	a.push_back(5);
	for(i = a.begin(); i!=a.end(); ++i){
		std::cout<<(*i)<<std::endl;
	}
	// begin() returns an iterator that "points" to the first element in the sequence
	// end() returns an iterator that "points" to one-past-the-last-element in the sequence
	// ++i used to advance from one element to the next
	return 0;
}
