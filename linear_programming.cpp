#include "linear_programming.hpp"


linear_programming::linear_programming(long double **A, long double *b, long double *c, unsigned c_size, unsigned b_size)
{
	this->c_size = c_size;
	this->c = new long double[c_size];
	this->b_size = b_size;
	this->t_lines = b_size + 1;
	this->t_columns = c_size + 1 + (2 * b_size);
	for (unsigned i = 0; i < c_size; i++)
		this->c[i] = c[i];
	this->tableau_t = make_tableau(A, b);
	this->basis_size = this->t_columns - this->b_size - 1;
	this->basis = new bool[this->basis_size]();
}

long double** linear_programming::make_tableau(long double **A, long double *b)
{
	for (unsigned i = 0; i < this->b_size; i++)
	{
		if (b[i] < 0)
		{
			this->t_columns += this->b_size;
			return auxiliary(A, b);
		}
	}
	return tableau(A, b);
}

void linear_programming::auxiliary_format(long double **t)
{
	for (unsigned i = 1; i < this->t_lines; i++)
	{
		if (t[i][this->t_columns - 1] < 0)
		{
			for (unsigned j = 0; j < 2 * this->b_size + this->c_size; j++)
			{
				if (t[i][j] != 0)
				{
					t[i][j] = -t[i][j];
					if (abs(t[i][j]) < THRESHOLD)
						t[i][j] = 0;
				}
			}
			t[i][this->t_columns - 1] = -t[i][this->t_columns - 1];
		}
		for (unsigned j = 0; j < this->t_columns; j++)
		{
			t[0][j] += -t[i][j];
			if (abs(t[0][j]) < THRESHOLD)
				t[0][j] = 0;
		}
	}
}

long double** linear_programming::auxiliary(long double **A, long double *b)
{
	long double **t;
	this->is_auxiliary = true;
	t = new long double*[this->t_lines]();
	for (unsigned i = 0, k1 = 0; i < this->t_lines; i++)
	{
		t[i] = new long double[this->t_columns]();
		if (i == 0)
		{
			for (unsigned j = 2 * this->b_size + this->c_size; j < this->t_columns - 1; j++)
				t[0][j] = 1;
		}
		else
		{
			for (unsigned j = this->b_size, k = 0; k < this->c_size; j++, k++)
			{
				t[i][j] = A[i - 1][k];
			}
			t[i][k1] = 1;
			t[i][k1 + this->b_size + this->c_size] = 1;
			t[i][k1 + 2 * this->b_size + this->c_size] = 1;
			t[i][this->t_columns - 1] = b[i - 1];
			k1++;
		}
	}
	auxiliary_format(t);
	return t;
}

long double** linear_programming::tableau(long double **A, long double *b)
{
	long double **t;
	this->is_auxiliary = false;
	t = new long double*[this->t_lines];
	for(unsigned i = 0, k1 = 0; i < this->t_lines; i++)
	{
		t[i] = new long double[this->t_columns]();
		if (i == 0)
		{
			for(unsigned j = this->b_size, k = 0; k < this->c_size; j++, k++)
				if (this->c[k] != 0)
					t[0][j] = -this->c[k];
		}
		else
		{
			for (unsigned j = this->b_size, k = 0; k < this->c_size; j++, k++)
				t[i][j] = A[i - 1][k];
			t[i][k1] = 1;
			t[i][k1 + this->b_size + this->c_size] = 1;
			t[i][this->t_columns - 1] = b[i - 1];
			k1++;
		}
	}
	return t;
}

inline void linear_programming::verify_threshold(unsigned i, unsigned j)
{
	if (abs(this->tableau_t[i][j]) < THRESHOLD)
		this->tableau_t[i][j] = 0;
}

inline void linear_programming::sum(int i, int j, long double constant)
{
	for (unsigned k = 0; k < this->t_columns; k++)
	{
		this->tableau_t[i][k] += constant * this->tableau_t[j][k];
		verify_threshold(i, k);
	}
}

bool linear_programming::verify_c(unsigned *column)
{
	for (unsigned i = this->b_size; i < this->t_columns - 1; i++)
	{
		if (this->tableau_t[0][i] < 0)
		{
			(*column) = i;
			return true;
		}
	}
	return false;
}

bool linear_programming::verify_column(unsigned i)
{
	bool unique = false;
	if (this->tableau_t[0][i] != 0)
		return false;
	for (unsigned j = 1; j < this->t_lines; j++)
	{
		if (this->tableau_t[j][i] != 0 && this->tableau_t[j][i] != 1)
			return false;
		if (this->tableau_t[j][i] == 1)
			unique = !unique;
	}
	return unique;
}

void linear_programming::find_basis()
{
	unsigned basis_counter = 0;
	for (unsigned i = this->b_size, j = 0; i < this->t_columns - 1; i++, j++)
	{
		if (verify_column(i))
		{
			this->basis[j] = true;
			basis_counter++;
			if (basis_counter == this->c_size)
				break;
		}
		else
			this->basis[j] = false;
	}
}

bool linear_programming::choose_pivot(unsigned column, unsigned *p_line)
{
	long double min = numeric_limits<double>::infinity();
	for (unsigned i = 1; i < this->t_lines; i++)
	{
		if (this->tableau_t[i][column] > 0 && (this->tableau_t[i][this->t_columns - 1] / this->tableau_t[i][column]) < min)
		{
			min = this->tableau_t[i][this->t_columns - 1] / this->tableau_t[i][column];
			(*p_line) = i;
		}
	}
	if (min == numeric_limits<double>::infinity())
	{
		this->unbound_column = column;
		return false;
	}
	return true;
}

void linear_programming::pivot_auxiliary()
{
	unsigned p_line = 0;
	for (unsigned i = 0, j = this->b_size; i < this->basis_size; i++, j++)
		if (this->basis[i] && choose_pivot(j, &p_line))
			sum(0, p_line, this->tableau_t[0][j]);
}

void linear_programming::simplex()
{
	long double pivot = 0;
	unsigned column = 0, p_line = 0;
	while (verify_c(&column))
	{
		if (choose_pivot(column, &p_line))
		{
			pivot = this->tableau_t[p_line][column];
			for (unsigned i = 0; i < this->t_columns; i++)
			{
				this->tableau_t[p_line][i] = this->tableau_t[p_line][i] / pivot;
				verify_threshold(p_line, i);
			}
			for (unsigned i = 0; i < this->t_lines; i++)
				if (i != p_line)
					sum(i, p_line, -this->tableau_t[i][column]);
		}
		else
		{
			this->result = "ilimitada";
			return;
		}
	}
	if (this->is_auxiliary && this->tableau_t[0][this->t_columns - 1] < 0)
	{
		this->result = "inviavel";
	}
	else if (!this->is_auxiliary)
	{
		this->result = "otima";
	}
}

void linear_programming::print_basis()
{
	for (unsigned i = 0; i < this->basis_size; i++)
		cout << this->basis[i] << " ";
	cout << "\n";
}

void linear_programming::solve_auxiliary()
{
	simplex();
	find_basis();
	if (this->result != "inviavel")
	{
		for (unsigned i = 0; i < this->b_size; i++)
			this->tableau_t[0][i] = 0;
		for (unsigned i = this->b_size, j = 0; i < this->b_size + this->c_size; i++, j++)
		{
			if (this->c[j] == 0)
				this->tableau_t[0][i] = 0;
			else
				this->tableau_t[0][i] = -c[j];
		}
		for (unsigned i = this->b_size + this->c_size; i < this->t_columns - this->b_size; i++)
			this->tableau_t[0][i] = 0;
		this->tableau_t[0][this->t_columns - 1] = 0;
		this->is_auxiliary = false;
		pivot_auxiliary();
		simplex();
	}
}

long double* linear_programming::find_solution_vector()
{
	unsigned p_line = 0;
	auto z = new long double[this->basis_size];
	find_basis();
	for (unsigned i = 0; i < this->basis_size; i++)
		if (this->basis[i]  && choose_pivot(this->b_size + i, &p_line))
		{
			z[i] = this->tableau_t[p_line][this->t_columns - 1];
			if (abs(z[i]) < THRESHOLD)
				z[i] = 0;
		}
		else
			z[i] = 0;
	return z;
}

long double* linear_programming::find_unbound_certificate()
{
	unsigned p_line = 0;
	auto z = new long double[this->t_columns - this->b_size - 1]();
	for (unsigned i = 0, j = this->b_size; i < this->t_columns - this->b_size - 1 && j < this->t_columns - 1; i++, j++)
	{
		if (this->basis[i] && choose_pivot(j, &p_line))
		{
			if (this->tableau_t[p_line][this->unbound_column] == 0)
				z[i] = 0;
			else
				z[i] = -this->tableau_t[p_line][this->unbound_column];
			if (z[i] < THRESHOLD)
				z[i] = 0;
		}
		else if (j == this->unbound_column)
			z[i] = 1;
	}
	return z;
}

void linear_programming::solve()
{
	if (this->is_auxiliary)
		solve_auxiliary();
	else
		simplex();
}
void linear_programming::print_solution()
{
	long double *z;
	cout << this->result << endl;
	if (this->result == "otima")
	{
		z = find_solution_vector();
		cout << this->tableau_t[0][this->t_columns - 1] << "\n";
		for (unsigned i = 0; i < this->c_size; i++)
			cout << z[i] << " ";
		cout << "\n";
		for (unsigned i = 0; i < this->b_size; i++)
			cout << this->tableau_t[0][i] << " ";
		cout << "\n";
		delete[] z;
	}
	else if (this->result == "inviavel")
	{
		for (unsigned i = 0; i < this->b_size; i++)
			cout << this->tableau_t[0][i] << " ";
		cout << "\n";
	}
	else if (this->result == "ilimitada")
	{
		z = find_solution_vector();
		for (unsigned i = 0; i < this->c_size; i++)
			cout << z[i] << " ";
		cout << "\n";
		delete[] z;
		z = find_unbound_certificate();
		for (unsigned i = 0; i < this->c_size; i++)
			cout << z[i] << " ";
		cout << "\n";
		delete[] z;
	}
}
void linear_programming::print_tableau()
{
	for (unsigned i = 0; i < this->t_lines; i++)
	{
		for (unsigned j = 0; j < this->t_columns; j++)
			cout << setw(5) << this->tableau_t[i][j] << " ";
		cout << "\n";
	}
}

linear_programming::~linear_programming()
{
	for(unsigned i = 0; i < this->t_lines; i++)
		delete[] this->tableau_t[i];
	delete[] this->tableau_t;
	delete[] this->basis;
	delete[] this->c;
}