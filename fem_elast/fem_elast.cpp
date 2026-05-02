#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

struct Node
{
	double x, y;
	bool isBoundary;
};

struct Segment
{
	double a, b;
};

struct FiniteElement
{
	std::vector<Node> nodes;
	std::vector<int> numNodes;
};

std::vector<Node> Grid(Segment& xSeg, Segment& ySeg, int& nx, int& ny)
{
	int gSize = (nx + 1) * (ny + 1);
	std::vector<Node> grid(gSize);
	double hx = (xSeg.b - xSeg.a) / nx;
	double hy = (ySeg.b - ySeg.a) / ny;

	for (int j = 0; j < ny + 1; ++j)
	{
		for (int i = 0; i < nx + 1; ++i)
		{
			int k = j * (nx + 1) + i;
			grid[k].x = xSeg.a + i * hx;
			grid[k].y = ySeg.a + j * hy;
			grid[k].isBoundary = i == 0 || i == nx || j == 0 || j == ny;
		}
	}
	return grid;
}

double dx(double& x, double& y, double(*f) (double, double), double eps = 1e-6)
{
	return (f(x + eps, y) - f(x, y)) / eps;
}

double dy(double& x, double& y, double(*f) (double, double), double eps = 1e-6)
{
	return (f(x, y + eps) - f(x, y)) / eps;
}

double dot(std::vector<double>& v, std::vector<double>& u)
{
	double s = 0;
	for (int i = 0; i < u.size(); ++i)
	{
		s += v[i] * u[i];
	}
	return s;
}

std::vector<double> grad(double& x, double& y, double(*f) (double, double))
{
	std::vector<double> v(2);
	v[0] = dx(x, y, f);
	v[1] = dy(x, y, f);
	return v;
}

int Delta(int& i, int& j)
{
	if (i == j)
		return 1;
	return 0;
}

std::vector<double> rhs1()
{
	return { 0, -10 };
}

int main()
{
	// Генерация сетки
	Segment xSeg{ 0, 10 };
	int nx = 100;
	Segment ySeg{ 0, 5 };
	int ny = 100;
	std::vector<Node> grid = Grid(xSeg, ySeg, nx, ny);

	// Константы
	double E = 21e10;
	double nu = 0.3;
	double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
	double mu = E / (2 * (1 + nu));

	// Создание конечных элементов
	std::vector<FiniteElement> finEls;
	for (int j = 0; j < ny; ++j)
	{
		for (int i = 0; i < nx; ++i)
		{
			std::vector<int> inds{ j * (nx + 1) + i, (j + 1) * (nx + 1) + i, (j + 1) * (nx + 1) + i + 1, j * (nx + 1) + i + 1 };
			FiniteElement feL, feU;
			feU.nodes = { grid[inds[0]], grid[inds[1]], grid[inds[2]] };
			feU.numNodes = { inds[0], inds[1], inds[2] };
			feL.nodes = { grid[inds[0]], grid[inds[2]], grid[inds[3]] };
			feL.numNodes = { inds[0], inds[2], inds[3] };
			finEls.push_back(feU);
			finEls.push_back(feL);
		}
	}

	// Матрицы жёсткости
	int gSize = (nx + 1) * (ny + 1);
	std::vector<std::vector<double>> K(2 * gSize, std::vector<double>(2 * gSize));
	std::vector<double> F(2 * (nx + 1) * (ny + 1));
	double hx = (xSeg.b - xSeg.a) / nx;
	double hy = (ySeg.b - ySeg.a) / ny;
	double area = hx * hy / 2;
	std::vector<std::vector<double>> C(3, std::vector<double>(3));
	C[0][0] = lambda + 2 * mu;
	C[1][1] = lambda + 2 * mu;
	C[0][1] = lambda;
	C[1][0] = lambda;
	C[2][2] = mu;
	std::vector<double> bU{ 0, -1 / hx, 1 / hx };
	std::vector<double> cU{ -1 / hy, 1 / hy, 0 };
	std::vector<double> bL{ -1 / hx, 0, 1 / hx };
	std::vector<double> cL{ 0, 1 / hy, -1 / hy };
	std::vector<double> f = rhs1();
	for (auto& elem : finEls)
	{
		// Формирование локальной матрицы жёсткости
		std::vector<std::vector<double>> B(3, std::vector<double>(6));
		// верхнетреугольный КЭ
		if (elem.numNodes[2] - elem.numNodes[1] == 1)
		{
			for (int i = 0; i < 3; ++i)
			{
				B[0][2 * i] = bU[i];
				B[1][2 * i + 1] = cU[i];
				B[2][2 * i] = cU[i] / sqrt(2);
				B[2][2 * i + 1] = bU[i] / sqrt(2);
			}
		}
		// нижнетреугольный КЭ
		else
		{
			for (int i = 0; i < 3; ++i)
			{
				B[0][2 * i] = bL[i];
				B[1][2 * i + 1] = cL[i];
				B[2][2 * i] = cL[i] / sqrt(2);
				B[2][2 * i + 1] = bL[i] / sqrt(2);
			}
		}
		std::vector<std::vector<double>> Ke(6, std::vector<double>(6));
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				double el = 0;
				for (int k = 0; k < 3; ++k)
				{
					for (int m = 0; m < 3; ++m)
					{
						el += B[k][i] * C[k][m] * B[m][j];
					}
				}
				Ke[i][j] = area * el;
			}
		}
		std::vector<double> Fe(6);
		for (int i = 0; i < 3; ++i)
		{
			Fe[2 * i] = f[0] * area / 3;
			Fe[2 * i + 1] = f[1] * area / 3;
		}

		// Внедрение Ke в K
		std::vector<int> globInd(6);
		for (int i = 0; i < 3; ++i)
		{
			globInd[2 * i] = 2 * elem.numNodes[i];
			globInd[2 * i + 1] = 2 * elem.numNodes[i] + 1;
		}
		for (int i = 0; i < 6; ++i)
		{
			int I = globInd[i];
			F[I] += Fe[i];
			for (int j = 0; j < 6; ++j)
			{
				int J = globInd[j];
				K[I][J] += Ke[i][j];
			}
		}
	}


	return 0;
}