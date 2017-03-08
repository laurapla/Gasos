#include <vector>

int main()
{
	std::vector<int> v;
	v.push_back(999);
	std::vector<int>::reverse_iterator r = v.begin();
	std::vector<int>::iterator i = r.base();	// will point to the last
												// element in the sequence
}
