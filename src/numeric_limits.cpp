/*! To determine numerical characteristics of the system.
 */

#include <iostream>

#include <limits>

// for INT_MIN etc
#include <climits>

// for 64 bit integer int64_t
#include <cstdint>

using namespace std;

int main()
{
	cout << numeric_limits<double>::epsilon() << " size " << sizeof(double) << endl;
	cout << numeric_limits<long double>::epsilon() << " size " << sizeof(long double) << endl;
	cout << INT_MIN << " " << INT_MAX << std::endl;
	cout << sizeof(int64_t);
	cout << endl;
	return 0;
}
