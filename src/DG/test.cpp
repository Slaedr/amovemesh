#include <iostream>
#include <vector>

using namespace std;

struct A
{
	int a;
	int b;
};

int main()
{
	vector<int> A;
	A.push_back(3); A.push_back(2); A.push_back(5);
	for(int i = 0; i < A.size(); i++)
		cout << A[i];
	cout << endl;
	A.erase(A.begin()+1);
	for(int i = 0; i < A.size(); i++)
		cout << A[i];
	cout << endl;

	cout << endl;
	return 0;
}
