
/* This programm performs Gaussian elemination with partial pivoting 
that solves given system of linear equations */




#include"system1.h"

int main()
{
	std::clock_t time = std::clock();

	try
	{
		int n;
		System A;
		std::cout << "Choose system input method (file/formula) ";
		std::string query;
		if ((std::cin >> query).fail()) throw std::logic_error ("Incorrect query");

		if (!((query == "file") || (query == "formula"))) throw std::logic_error ("Incorrect query");

		if (query == "formula") 
		{
			std::cout << "Insert system dimension ";
			if ((std::cin >> n).fail()) throw std::logic_error ("Incorrect query");
			A.Formula(n);
		}

		else n = A.File();

		//Set of variables
		double* X = new double [n];
	
		//Source system is needed to calculate the residual
		System B = A;

		A.Gauss(X);

		//The number of solution set elements that will be printed on the screen
		int m;

		std::cout << "Insert m ";
		std::cin >> m;
		B.Print(X, m);
		delete [] X;
	}
	catch (std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}

	time = (std::clock() - time)/(double)CLOCKS_PER_SEC;
	std::cout << std::fixed;
	std::cout.precision(4);
	std::cout << "Elapsed time: " << time << " seconds\n";

	return 0;
}
