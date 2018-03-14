#include<iostream>
#include<string>
#include<fstream>
#include<stdexcept>
#include<cmath>
#include<ctime>
#include<iomanip>
#include <cstdlib>

//Measurement error
#define EPSILON 1e-8

//Random function to set the matrix
double f (int i, int j);
double f (int i, int j)
{
	double a;
	a = pow(-1, i*j)*(i + 1)/(j + 1);
	return a;
}

class System
{
	private:
		//Coefficient matrix
		double** A;

		//Column vector of constant terms
		double* B;

		//Permutation array
		int* I;

		//System dimension
		int n;
	public:
		//Constructor
		System ()
		{
			n = 0;
			A = NULL;
			B = NULL;
			I = NULL;
		}

		//Copy constructor
		System (const System &obj)
		{
			n = obj.n;
			A = new double* [n];			
			B = new double [n];
			I = new int [n];
			for (int i = 0; i < n; i++)
			{
				A[i] = new double [n];
				B[i] = obj.B[i];
				I[i] = obj.I[i];
				for (int j = 0; j < n; j++) A[i][j] = obj.A[i][j];
			}
		}

		//Sets system from file "input.txt"
		int File ()
		{	
			std::fstream file;
			file.open ("input.txt");
			if (!(file.is_open()))
			{
				throw std::logic_error ("File error");
			}
			if ((file >> n).fail()) throw std::logic_error ("Incorrect source matrix");
			A = new double* [n];
			B = new double [n];
			I = new int [n];
			for (int i = 0; i < n; i++)
			{
				I[i] = i;
				A[i] = new double [n];
				for (int j = 0; j < n; j++)
				{
					if ((file >> A[i][j]).fail()) throw std::logic_error ("Incorrect source matrix");
				}
				if ((file >> B[i]).fail()) throw std::logic_error ("Incorrect source matrix");
			} 
			file.close();
			return n;
		}

		//Sets system using function f
		void Formula (int m)
		{
			n = m;
			A = new double* [n];
			B = new double [n];
			I = new int [n];
			for (int i = 0; i < n; i++)
			{
				A[i] = new double [n];
				I[i] = i;
				double s = 0;
				for (int j = 0; j < n; j++)
				{
					A[i][j] = f(i, j);
					if (j%2 == 1) s = s + A[i][j];
				}
				B[i] = s;
			}
		}

		//Destructor
		~System ()
		{
			for (int i = 0; i < n; i++) delete A[i];
			delete [] A;
			delete [] B;
			delete [] I;
		}

		//Prints the system on the screen
		void Show()
		{
			if (!(n > 0)) std::cout << "Matrix is empty" << std::endl;
			else
			{
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						if (j < n - 1)
						{
							if (A[i][j + 1] >= 0) std::cout << A[i][j] << "x" << I[j] + 1 << "+";
							else std::cout << A[i][j] << "x" << I[j] + 1;
						}
						else std::cout << A[i][j] << "x" << I[j] + 1 << " = " << B[i] << std::endl;
					}
				}
				std::cout << "______________________________________\n\n";
			}	
		}

		//Solves the system using Gaussian elemination with partial pivoting and assigns the solution set to array X
		void Gauss (double X[])
		{
			if (!(n > 0)) std::cout << "Matrix is empty" << std::endl;
			else
			{
 				//The nubmer of null rows
				int h = 0;

				//array Q is used later to mark a null row (Q[i] = -1 if the row number i is null)
				int* Q = new int [n];
				for (int i = 0; i < n; i++)
				{
					//The column of leading element
					int l = i - h;

					//Finds the greatest in absolute value element in row number i
					double pivot = 0;
					for (int j = i - h; j < n; j++)
					{
						if (fabs(A[i][j]) > pivot)
						{
							pivot = fabs(A[i][j]);
							l = j;
						}
					}

					//If the greatest by absolute value element equals zero
					//then It is a null row. Mark it assigning value -1 to Q[i]
					if (fabs(pivot) < EPSILON*n)						
					{									
						Q[i] = -1;							
						h++;
						
						//If corresponding constant term is not equal to zero then the systen is inconsistent
						if (fabs(B[i]) > EPSILON*n)
						{								
							delete [] Q;						
							delete [] X;						
							throw std::logic_error ("System is inconsisted");
						}
					}

					//Else (if the greatest by absolute value element is not equal to zero)
					//start Gaussian elemination with partial pivoting
					else							
					{							
						Q[i] = h;					
						if (l != (i - h))				
						{		
							//Swaps columns if needed
							for (int j = 0; j < n; j++)
							{
								double m = A[j][l];
								A[j][l] = A[j][i - h];
								A[j][i - h] = m;
							}

							//Remembers the swap using the permutation array
							int g = I[i - h];
							I[i - h] = I[l];
							I[l] = g;
						}

						//Divides the row by the value of its leading element
						double leader = A[i][i - h];
						for (int j = i - h; j < n; j++) A[i][j] = A[i][j]/leader;
						B[i] = B[i]/leader;

						//Eleminates the lower rows
						for (int p = i + 1; p < n; p++)
						{
							double d = A[p][i - h];
							B[p] = B[p] - B[i]*d;
							for (int j = i - h; j < n; j++) A[p][j] = A[p][j] - A[i][j]*d;
						}
					}
				}

				//If at least one row is null and the corresponding constant terms for all null rows equal zero
				//then the system has infinitely many solutions, i.e. it has free variables.
				//We assign value 1 to all free variables
				for (int j = 0; j < h; j++) X[I[n - 1 - j]] = 1;


				//Back substitution
				for (int i = n - 1; i >= 0; i--)
				{
					if (Q[i] > -1)
					{ 
						double s = B[i];
						for (int k = i - Q[i] + 1; k < n; k++) s = s - A[i][k]*X[I[k]];
						X[I[i - Q[i]]] = s;
					}
				}
				delete [] Q;
			}
		}

		//Prints residual and the first m elements of the solution set on the screen 
		void Print (double X[], int m)
		{
			if (m > n)
			{
				for (int i = 0; i < n; i++) std::cout << "x" << i + 1 << "=" << X[i] << std::endl;
			}
			else
			{
				for (int i = 0; i < m; i++) std::cout << "x" << i + 1 << "=" << X[i] << std::endl;
			}
			double res = 0;
			for (int i = 0; i < n; i++)
			{
				double s = 0;
				for (int j = 0; j < n; j++) s = s + A[i][j]*X[j];
				res = res + (s - B[i])*(s - B[i]);
			}
			res = sqrt(res);
			std::cout << "Residual = " << res << std::endl;
		}
};
