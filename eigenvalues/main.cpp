
/*This programm finds eigenvalues (if they are real) of given matrix using QR decomposition*/




#include"qr.h"

int main()
{
	std::clock_t time = clock();

	try {	
		int n;
		Matrix A;
		std::cout << "Choose system input method (file/formula) ";
		std::string str;
		std::cin >> str;
		if (!((str == "file") || (str == "formula"))) throw std::logic_error ("Incorrect query");

		if (str == "formula")
		{
			std::cout << "Insert matrix dimension ";
			std::cin >> n;
			A.Formula(n);
		}
		else n = A.File();

		double* L = new double [n];
		Matrix B = A;
		double Im;
		Im = A.QR(L);
		int m;
		std::cout << "Insert m ";
		std::cin >> m;
		A.Print(m, B, L, Im);
		
	}
	catch(std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	time = (std::clock() - time)/(double)CLOCKS_PER_SEC;
	std::cout << "Elapsed time: " << time << " seconds\n";

	return 0;
}
