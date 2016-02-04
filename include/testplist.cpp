#include "adatastructures.hpp"

using namespace std;

int main()
{
	PList<double> list;
	list.push_back(3.14);
	list.push_back(2.478);
	cout << list[0] << endl;
	list[0] = 5.67;
	list.push_back(9.99);
	list.print();
	cout << "*\n";
	
	list.delete_element(1);
	list.print();

	cout << endl;
	return 0;
}
