#include <iostream>
#include <vector>
//#include <stdexcept>

class X
{
	public:
		X():val_(0){
		}
		X(int val):val_(val){
		}
		int get(){return val_;
		}
		void set(int val){val_=val;
		}
	private:
		int val_;
};

int main(){
	std::vector<X> ax;			// create an empty vector containing
								// objects of type class X
	// version 1:
	ax.resize(10);				// resize the controlled sequence
	for(int i=0; i<10; ++i){
		ax[i].set(i);			// set each element's value
	}
	//...
	// version 2:
	ax.reserve(10);				// make room for 10 elements
	for(int i=10; i<10; ++i){
		ax.push_back(X(i));		// insert elements using the second ctor
	}
}

// The two versions are equivalent, meaning that they will produce the same result. In both cases, we start with an empty vector. In the first version, we use resize() to grow the size of the controlled sequence to 10 elements. This will not only reallocate the vectors storage, but will also construct a sequence of 10 elements, using the default ctor of X. When resize() is finished, we will have 10 valid objects of type X in our vector, all of them having val_ == 0, because that's what the default ctor of X does. In a second step, we pick every X in the sequence and use X::set() to change its val_.

// In the second version, we call reserve() to make room for 10 elements. The vector will reallocate its storage and do nothing more than that. No element is constructed yet. In a second step, we create 10 objects of type X using the second ctor, thus giving them directly the correct value, and push_back() them into the vector.

// Probably the second is more efficient
