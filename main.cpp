#include "linear_programming.hpp"

int main()
{
	unsigned int i = 0, j = 0;
	long double **A, *b, *c;
	string s, token;
	getline(cin, s);
	unsigned b_size = stoi(s.substr(0, s.find(' ')));
	unsigned c_size = stoi(s.substr(s.find(' ')));
	b = new long double[b_size]();
	c = new long double[c_size]();
	A = new long double*[b_size];
	getline(cin, s);
	istringstream ss(s);
	while(getline(ss, token, ' '))
	{
		c[i] = stold(token);
		i++;
	}
	for(i = 0; i < b_size; i++)
	{
		getline(cin, s);
		istringstream ss(s);
		A[i] = new long double[c_size]();
		while(getline(ss, token, ' '))
		{
			if (j == c_size)
			{
				b[i] = stold(token);
				break;
			}
			A[i][j] = stold(token);
			j++;
		}
		j = 0;
	}
	auto lp = linear_programming(A, b, c, c_size, b_size);
	delete[] c;
	delete[] b;
	for (unsigned i = 0; i < b_size; i++)
	{
		delete[] A[i];
	}
	delete[] A;
	lp.solve();
	lp.print_solution();
	return 0;
}