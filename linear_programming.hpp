#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <string>
#define THRESHOLD 0.00000000000001
using namespace std;

class linear_programming
{
	private:
		long double **tableau_t;
		long double *c;
		bool *basis;
		bool is_auxiliary;
		string result;
		unsigned c_size;
		unsigned b_size;
		unsigned basis_size;
		unsigned unbound_column;
		unsigned t_lines;
		unsigned t_columns;
		long double **make_tableau(long double **A, long double *b);
		long double **tableau(long double **A, long double *b);
		void auxiliary_format(long double **t);
		long double **auxiliary(long double **A, long double *b);
		void verify_threshold(unsigned i, unsigned j);
		void sum(int i, int j, long double constant = 1);
		void pivot_auxiliary();
		bool choose_pivot(unsigned column, unsigned *p_line);
		bool verify_column(unsigned i);
		bool verify_c(unsigned *column);
		long double *find_solution_vector();
		long double *find_unbound_certificate();
		void find_basis();
		void simplex();
		void solve_auxiliary();
	
	public:
		linear_programming(long double **A, long double *b, long double *c, unsigned c_size, unsigned b_size);
		void solve();
		void print_solution();
		void print_tableau();
		void print_basis();
		~linear_programming();
};