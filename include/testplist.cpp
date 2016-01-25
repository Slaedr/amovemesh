#include "adatastructures.hpp"

using namespace std;

int main()
{
	PList<double> list;
	list.push_back(3.14);
	list.push_back(2.478);
	cout << list[0] << endl;
	list[0] = 5.67;
	cout << list[0] << " " << list[1] << endl;
	list.print();

	cout << endl;
	return 0;
}
