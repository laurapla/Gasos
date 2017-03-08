#include <iostream>
#include <vector>
//#include <stdexcept>


int main(){
	std::vector<char> array;
	int i = 999;		// some integer value
	array.reserve(10);	// make room for 10 elements
	array.push_back(i);
	std::cout<<array.capacity()<<std::endl;
	std::cout<<array.size()<<std::endl;
}
