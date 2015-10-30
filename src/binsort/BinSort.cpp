#include<stdio.h>
#include<math.h>
#include<iostream>
#include <fstream>
#include <string>
using namespace std;

int main(){
	string filename;
	filename = "bounpoints.dat";
	ifstream fin("bounpoints.dat");
	if (!fin)
	{
		cout << "error" << endl;
		return 0;
	}
	int npoints = 499;
	double points[499][4];
	for (int i = 0; i < 499; i++)
	{
		for (int j = 0; j < 2; j++)
			fin >> points[i][j];
	}
	double xmin = points[0][0], xmax = points[0][0], ymin = points[0][1], ymax = points[0][1];
	int rownum = 2; int colnum = 24;
	for (int i = 1; i < 499; i++)
	{
		if (points[i][0] > xmax) xmax = points[i][0];
		if (points[i][0] < xmin) xmin = points[i][0];
		if (points[i][1] > ymax) ymax = points[i][1];
		if (points[i][1] < ymin) ymin = points[i][1];
	}
	xmax = ceil(xmax) + 1; xmin = floor(xmin) - 1;
	ymax = ceil(ymax) + 1; ymin = floor(ymin) - 1;
	for (int i = 0; i < npoints; i++)
	{
		int ysteps = floor((points[i][1] - ymin)*rownum / (ymax - ymin));
		points[i][2] = ysteps*colnum + (colnum - 1 + int(pow(-1, ysteps + 1))*(colnum - 1 - 2 * floor((points[i][0] - xmin)*colnum / (xmax - xmin)))) / 2;
	}
	int pointSeq[499];
	int index = 0;
	for (int i = 0; i < npoints; i++)
	{
		points[i][3] = 0;
	}
	for (int k =12; k < 24; k++)
	{
		printf("key=%d\n", k);
		for (int i = 0; i < npoints; i++)
		{
			if ((points[i][2] == k) && (points[i][3] == 0))
			{
				pointSeq[index] = i;
				points[i][3] = 1;
				index++;
				printf("p%d\n", i);
			}
		}
		printf("\n");
	}
	ofstream fin2("sortingData.txt");
	for (int i = 0; i < 499; i++)
	{
		for (int j = 0; j < 2; j++)
			fin2 << points[pointSeq[i]][j] << " ";
		fin2 << '\n';
	}
	
	return 0;
}
