#include <vector>

std::vector<int> v;
v.push_back(999);
// fill up the vector

int main()
{
	// following statements are equivalent:
	int i = v.front();
	int i = v[0];
	int i = v.at(0);
	int i = *(v.begin());
	// following statements are equivalent:
	int j = v.back();
	int j = v[v.size()-1];
	int j = v.at(v.size()-1);
	int j = *(v.end()-1);
}
