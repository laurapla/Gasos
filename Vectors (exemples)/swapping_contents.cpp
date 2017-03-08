#include <vector>

std::vector<int> v;

int main()
{
	v.clear();
	v.swap(std::vector<int>(v));
}
