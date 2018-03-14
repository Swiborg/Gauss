#include<iostream>
#include<string>
#include<fstream>
#include<stdexcept>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include <algorithm>

//Measurement error
#define EPSILON 1e-8

double f (int i, int j);
double f (int i, int j)
{
	double a;
	a = pow(-1, i*j)*(i + 1)/(j + 1);
	return a;
}

class Matrix
{
	private:
		//Matrix
		double** A;
		
		//Size of matrix
		int n;
	public:

		Matrix ()
		{
			n = 0;
			A = NULL;
		}

		Matrix (int m)
		{
			n = m;
			A = new double* [n];
			for (int i = 0; i < n; i++)
			{
				A[i] = new double [n];
			}
		}

		Matrix (const Matrix &obj)
		{
			n = obj.n;
			A = new double* [n];
			for (int i = 0; i < n; i++)
			{
				A[i] = new double [n];
				for (int j = 0; j < n; j++)
				{
					A[i][j] = obj.A[i][j];
				}
			}
		}

		void E (int m){
			n = m;
			A = new double* [n];
			for(int i = 0; i < n; i++)
			{
				A[i] = new double [n];
			}
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					if (i == j) A[i][j] = 1;
					if (i != j) A[i][j] = 0;
				}
			}
		}
		
		//Sets system from file "input.txt"
		int  File ()
		{
			std::fstream file;
			file.open ("input.txt");
			if (!(file.is_open())) throw std::logic_error ("File error");
			if ((file >> n).fail()) throw std::logic_error ("Incorrect source matrix");
			A = new double* [n];
			for (int i = 0; i < n; i++)
			{
				A[i] = new double [n];
				for(int j = 0; j < n; j++)
				{
					if ((file >> A[i][j]).fail()) throw std::logic_error ("Incorrect source matrix");
				}
			} 
			file.close();
			return n;
		}
		
		//Sets system using function f
		void Formula (int m)
		{
			n = m;
			A = new double* [n];
			for (int i = 0; i < n; i++)
			{
				A[i] = new double [n];
				for (int j = 0; j < n; j++)
				{
					A[i][j] = f(i,j);
				}
			}
		}

		~Matrix()
		{
			for (int i = 0; i < n; i++) delete A[i];
			delete [] A;
		}

		void Show()
		{
			if (!(n > 0)) std::cout << "Matrix is empty\n";
			else
			{
				for(int i = 0; i < n; i++)
				{
					for(int j = 0;j < n; j++)
					{
						std::cout << A[i][j] << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "______________________________________\n\n";
			}	
		}

		void Rotation()
		{
			for (int j = 0; j < n - 2; j++)
			{
				for (int i = j + 2; i < n; i++)
				{
						if ((A[i][j] > 0) || (A[i][j] < 0))
						{	
							double x = A[j + 1][j] / sqrt(A[j + 1][j]*A[j + 1][j] + A[i][j]*A[i][j]);
							double y = A[i][j] / sqrt(A[j + 1][j]*A[j + 1][j] + A[i][j]*A[i][j]);
							for (int k = j; k < n; k++)
							{
								double a = A[j + 1][k], b = A[i][k];
								A[j + 1][k] = x*a + y*b;
								A[i][k] = -y*a + x*b;
							}
							for (int k = 0; k < n; k++)
							{
								double a = A[k][j + 1], b = A[k][i];
								A[k][j + 1] = x*a + y*b;
								A[k][i] = -y*a + x*b;
							}
						}
				}
			}
		}

		void Reflexion (Matrix *Q, int Iteration)
		{
			Q->E(Iteration + 1);
			for (int j = 0; j < Iteration; j++)
			{
				if ((A[j + 1][j] < 0) || (A[j + 1][j] > 0))
				{
					double x, y, z, r;
					r = sqrt(A[j][j]*A[j][j] + A[j + 1][j]*A[j + 1][j]);
					x = 1 - 2*(A[j][j] - r)*(A[j][j] - r)/((A[j][j] - r)*(A[j][j] - r) + A[j + 1][j]*A[j + 1][j]);
					y = -2*(A[j][j] - r)*A[j + 1][j]/((A[j][j] - r)*(A[j][j] - r) + A[j + 1][j]*A[j + 1][j]);
					z = 1 - 2*A[j + 1][j]*A[j + 1][j]/((A[j][j] - r)*(A[j][j] - r) + A[j + 1][j]*A[j + 1][j]);
					for (int k = 0; k <= Iteration; k++)
					{
						double a = Q->A[k][j], b = Q->A[k][j + 1];
						Q->A[k][j] = x*a + y*b;
						Q->A[k][j + 1] = y*a + z*b;
					}
					for(int i = j; i <= Iteration; i++)
					{
						double a = A[j][i], b = A[j + 1][i];
						A[j][i] = x*a + y*b;
						A[j + 1][i] = y*a + z*b;
					}
				}
			}
		}

		double QR (double L[])
		{
			double Im = 0;
			if (!(n > 0)) throw std::logic_error ("Matrix is empty");
			if (n == 1) throw std::logic_error ("λ1=" + std::to_string(A[0]));
			else
			{
				Matrix Q(n);
				Rotation();
				Matrix T(n);
				for (int Iteration = n - 1; Iteration > 1; Iteration--)
				{
					double v = 0;
					while (fabs(A[Iteration][Iteration - 1]) > EPSILON)
					{
						double w = A[n - 1][n - 1];
						for (int i = 0; i <= Iteration; i++)
						{
							A[i][i] = A[i][i] - w;
							for (int j = 0; j <= Iteration; j++)
							{
								if (fabs(A[i][j]) > v) v = fabs(A[i][j]);
							}
						}
						Reflexion(&Q, Iteration);
						for (int i = 0; i <= Iteration; i++)
						{
							for (int j = 0; j <= Iteration; j++)
							{
								double s = 0;
								for (int k = 0; k <= Iteration; k++) s = s + A[i][k]*Q.A[k][j];
								T.A[i][j] = s;
							}
						}
						for (int i = 0; i <= Iteration; i++)
						{
							for (int j = 0; j <= Iteration; j++)
							{
								if (i == j) A[i][j] = T.A[i][j] + w;
								else A[i][j] = T.A[i][j];
							}
						}
					}
					L[Iteration] = A[Iteration][Iteration];
				}	

				double D;
				D = (A[0][0] + A[1][1])*(A[0][0] + A[1][1]) - 4*(A[0][0]*A[1][1] - A[0][1]*A[1][0]);
				if (D < 0)
				{
					L[1] = (A[0][0] + A[1][1])/2;
					L[0] = (A[0][0] + A[1][1])/2;
					Im = sqrt(-D)/2;
				}
				else
				{
					L[1] = (A[0][0] + A[1][1] + sqrt(D))/2;
					L[0] = (A[0][0] + A[1][1] - sqrt(D))/2;
					Im = -1;
				}
			}
			return Im;
		}

		void Print (int m, Matrix B, double L[], double Im)
		{
			double a = 0, b = 0, c = 0, d = 0;
			if (Im < 0)
			{
				std::cout << "λ" << 1 << "=" << L[0] << std::endl;
				std::cout << "λ" << 2 << "=" << L[1] << std::endl;
			}
			else
			{
				std::cout << "λ" << 1 << "=" << L[0] << "-" << Im << "i" << std::endl;
				std::cout << "λ" << 2 << "=" << L[1] << "+" << Im << "i" << std::endl;
			}
			for (int i = 2; i < std::min(n, m); i++) std::cout << "λ" << i + 1 << "=" << L[i] << std::endl;
			for (int i = 0; i < n; i++)
			{
				a = a + A[i][i];
				b = b + B.A[i][i];
				for (int j = 0; j < n;j++)
				{
					c = c + A[i][j]*A[i][j];
					d = d + B.A[i][j]*B.A[i][j];
				}
			}
			a = fabs(a - b);
			b = fabs(sqrt(c) - sqrt(d));
			std::cout << "Trace residual: " << a << std::endl;
			std::cout << "Norm residual: " << b << std::endl;
		}
};
