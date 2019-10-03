#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;
string infile = "AMD.csv";
int n;
vector <vector <int> > f;	
int main(int agrv, char ** argc)
{	
	infile = argc[1];
	ifstream fi(infile.c_str());
	string s;
	int x;
	while (getline(fi,s))
	{
		// cerr << s << endl;
		// cerr << s.length() << endl;		
		for (int i = 0; i < s.length(); ++i)
		{
			if (s[i] == ',') cout << " ";
			else cout << s[i];
		}		
		cout << endl;
		++n;
	}
	fi.close();	
	return 0;
}