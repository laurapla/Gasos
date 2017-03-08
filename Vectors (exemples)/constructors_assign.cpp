#include <iostream>
#include <vector>
#include <list>

int main()
{
	typedef std::vector<std::string> str_vec_t;
	str_vec_t v1;						// create an empty vector
	str_vec_t v2(10);					// 10 copies of empty strings
	str_vec_t v3(10, "hello");			// 10 copies of the string
										// "hello"
	str_vec_t v4(v3);					//copy ctor
	
		std::list<std::string> sl;		// create a list of strings
										// and populate it
		sl.push_back("cat");
		sl.push_back("dog");
		sl.push_back("mouse");
		
	str_vec_t v5(sl.begin(), sl.end());	// a copy of the range in
										// another container
										// (here, a list)
	
	v1 = v5;							// will copy all elements
										// from v5 to v1
	
	// THE assign() FUNCTION: the old elements of the vector are discarded and the size of the vector is set to the number of elements assigned
	v1.assign(sl.begin(), sl.end());	// copies the list into
										// the vector
	v1.assign(3, "hello");				// initializes the vector
										// with 3 strings "hello"
	
}
