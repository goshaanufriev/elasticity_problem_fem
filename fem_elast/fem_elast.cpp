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
	int type{}; // 0 - L, 1 - U
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

double cubeNorm(std::vector<std::vector<double>>& A, const int& n)
{
	double norm = 0;
	for (int i = 0; i < n; ++i)
	{
		double s = 0;
		for (int j = 0; j < n; ++j)
		{
			s += std::fabs(A[i][j]);
		}
		if (s > norm)
			norm = s;
	}
	return norm;
}

double cubeNorm(std::vector<double>& b, const int& n)
{
	double norm = 0;
	for (int i = 0; i < n; ++i)
	{
		if (std::fabs(b[i]) > norm)
			norm = std::fabs(b[i]);
	}
	return norm;
}

std::vector<double> multMat(const std::vector<std::vector<double>>& A, const std::vector<double>& B)
{
	int n = B.size();
	std::vector<double> y(n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			y[i] += A[i][j] * B[j];
		}
	}
	return y;
}

std::vector<double> addVec(const std::vector<double>& x, const std::vector<double>& y)
{
	int n = x.size();
	std::vector <double> z(n);
	for (int i = 0; i < n; ++i)
	{
		z[i] = x[i] + y[i];
	}
	return z;
}

std::vector<double> subVec(const std::vector<double>& x, const std::vector<double>& y)
{
	int n = x.size();
	std::vector <double> z(n);
	for (int i = 0; i < n; ++i)
	{
		z[i] = x[i] - y[i];
	}
	return z;
}

//std::vector<double> Jacobi(std::vector<std::vector<double>>& A, std::vector<double>& b, const int& n)
//{
//	std::vector<std::vector<double>> C(n, std::vector<double>(n));
//	std::vector<double> y(n);
//	
//	for (int i = 0; i < n; ++i)
//	{
//		for (int j = 0; j < n; ++j)
//		{
//			if (i == j)
//			{
//				C[i][j] = 0;
//				y[i] = b[i] / A[i][j];
//			}
//			else
//			{
//				C[i][j] = -A[i][j] / A[i][i];
//			}
//		}
//	}
//
//	std::vector<double> x0(n);
//	for (int i = 0; i < n; ++i)
//	{
//		x0[i] = 0;
//	}
//
//	double eps = 1e-6;
//	double normC = cubeNorm(C, n);
//
//	std::vector<double> cx = multMat(C, x0, n);
//	std::vector<double> x = addVec(cx, y, n);
//	std::vector<double> err = subVec(x, x0, n);
//	double errNorm = cubeNorm(err, n);
//	std::cout << errNorm << "\n";
//	while (errNorm > eps)
//	{
//		x0 = x;
//		cx = multMat(C, x0, n);
//		x = addVec(cx, y, n);
//		err = subVec(x, x0, n);
//		errNorm = cubeNorm(err, n);
//		std::cout << errNorm << "\n";
//	}
//	return x;
//}

std::vector<double> kVec(std::vector<double>& v, double& k)
{
	std::vector<double> u(v.size());
	for (int i = 0; i < u.size(); ++i)
	{
		u[i] = k * v[i];
	}
	return u;
}

//std::vector<double> conjGrad(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x0)
//{
//	int n = b.size();
//	double alpha = 0;
//	std::vector<double> ax = multMat(A, x0);
//	std::vector<double> r = subVec(ax, b);
//	std::vector<double> ar = multMat(A, r);
//	double beta = -dot(r, r) / dot(ar, r);
//	std::vector<double> delta = kVec(r, beta);
//	std::vector<double> x = addVec(x0, delta);
//	std::vector<double> adx, dalpha, rbeta;
//	for (int i = 1; i < n; ++i)
//	{
//		ax = multMat(A, x);
//		r = subVec(ax, b);
//		ar = multMat(A, r);
//		adx = multMat(A, delta);
//		// dots
//		double rr = dot(r, r);
//		double ardx = dot(ar, delta);
//		double rdx = dot(r, delta);
//		double arr = dot(ar, r);
//		double adxdx = dot(adx, delta);
//
//		double denom = adxdx * arr - ardx * ardx;
//		double alpha = (rr * ardx - rdx * arr) / denom;
//		double beta = (rdx * ardx - rr * adxdx) / denom;
//		dalpha = kVec(delta, alpha);
//		rbeta = kVec(r, beta);
//		delta = addVec(dalpha, rbeta);
//		x = addVec(x, delta);
//	}
//	return x;
//}

std::vector<double> conjGrad(
	const std::vector<std::vector<double>>& A,
	const std::vector<double>& b,
	std::vector<double> x)
{
	int n = b.size();

	std::vector<double> r = subVec(b, multMat(A, x));
	std::vector<double> p = r;

	double rr = dot(r, r);
	double eps = 1e-14;

	for (int k = 0; k < n; ++k)
	{
		std::vector<double> Ap = multMat(A, p);

		double alpha = rr / dot(p, Ap);

		for (int i = 0; i < n; ++i)
			x[i] += alpha * p[i];

		for (int i = 0; i < n; ++i)
			r[i] -= alpha * Ap[i];

		double rr_new = dot(r, r);

		if (sqrt(rr_new) < eps)
			break;

		double beta = rr_new / rr;

		for (int i = 0; i < n; ++i)
			p[i] = r[i] + beta * p[i];

		rr = rr_new;
	}

	return x;
}

std::vector<double> rhs1()
{
	return { 0, -10 };
}

struct Coeffs
{
	double a, b, c;
};

std::vector<double> uEx(double& x, double& y, Coeffs& uxCoeff, Coeffs& uyCoeff)
{
	double a1 = uxCoeff.a, a2 = uxCoeff.b, a3 = uxCoeff.c;
	double b1 = uyCoeff.a, b2 = uyCoeff.b, b3 = uyCoeff.c;
	double ux = a1 * x * x + a2 * x * y + a3 * y * y;
	double uy = b1 * x * x + b2 * x * y + b3 * y * y;
	return { ux, uy };
}

int main()
{
	// Генерация сетки
	Segment xSeg{ 0, 10 };
	int nx = 10;
	Segment ySeg{ 0, 5 };
	int ny = 10;
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
			feU.type = 1;
			feL.nodes = { grid[inds[0]], grid[inds[2]], grid[inds[3]] };
			feL.numNodes = { inds[0], inds[2], inds[3] };
			feL.type = 0;
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
		if (elem.type == 1)
		{
			for (int i = 0; i < 3; ++i)
			{
				B[0][2 * i] = bU[i];
				B[1][2 * i + 1] = cU[i];
				B[2][2 * i] = cU[i];
				B[2][2 * i + 1] = bU[i];
			}
		}
		// нижнетреугольный КЭ
		else
		{
			for (int i = 0; i < 3; ++i)
			{
				B[0][2 * i] = bL[i];
				B[1][2 * i + 1] = cL[i];
				B[2][2 * i] = cL[i];
				B[2][2 * i + 1] = bL[i];
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

	// Учёт ГУ 1 рода
	// y = 0: (a1*x^2, b1*x^2)
	// x = lx: (a1*lx^2+a2*lx*y+a3*y^2, b1*lx^2+b2*lx*y+b3*y^2)
	// y = ly: (a1*x^2+a2*x*ly+a3*ly^2, b1*x^2+b2*x*ly+b3*ly^2)
	// x = 0: (ux=a3*y^2, uy=b3*y^2)

	double a1 = 0;
	double a3 = 1;
	double b1 = 1;
	double b3 = 0;
	double f1 = 0;
	double f2 = -1;

	double a2 = -(2 * mu * b1 + (2 * lambda + 4 * mu) * b3 + f2) / (lambda + mu);
	double b2 = -(2 * mu * a3 + (2 * lambda + 4 * mu) * a1 + f1) / (lambda + mu);

	double eps = 1e-10;
	for (int i = 0; i < gSize; ++i)
	{
		double x = grid[i].x;
		double y = grid[i].y;

		int I = 2 * i;
		int J = 2 * i + 1;

		double ux, uy;
		bool flag = false;
		if (std::fabs(y - ySeg.a) < eps)
		{
			ux = a1 * x * x;
			uy = b1 * x * x;
			flag = true;
		}
		else if (std::fabs(x - xSeg.b) < eps)
		{
			ux = a1 * xSeg.b * xSeg.b + a2 * xSeg.b * y + a3 * y * y;
			uy = b1 * xSeg.b * xSeg.b + b2 * xSeg.b * y + b3 * y * y;
			flag = true;
		}
		else if (std::fabs(y - ySeg.b) < eps)
		{
			ux = a1 * x * x + a2 * ySeg.b * x + a3 * ySeg.b * ySeg.b;
			uy = b1 * x * x + b2 * ySeg.b * x + b3 * ySeg.b * ySeg.b;
			flag = true;
		}
		else if (std::fabs(x - xSeg.a) < eps)
		{
			ux = a3 * y * y;
			uy = b3 * y * y;
			flag = true;
		}

		if (flag)
		{
			// ux
			for (int j = 0; j < 2 * gSize; ++j)
			{
				F[j] -= K[j][I] * ux;
				K[I][j] = 0;
				K[j][I] = 0;
			}
			K[I][I] = 1;
			F[I] = ux;

			// uy
			for (int j = 0; j < 2 * gSize; ++j)
			{
				F[j] -= K[j][J] * uy;
				K[J][j] = 0;
				K[j][J] = 0;
			}
			K[J][J] = 1;
			F[J] = uy;
		}
	}

	/*for (int i = 0; i < 2 * gSize; ++i)
	{
		for (int j = 0; j < 2 * gSize; ++j)
		{
			std::cout << K[i][j] << "\t";
		}
		std::cout << "\n";
	}*/

	int sysSize = 2 * gSize;
	std::vector<double> u0(sysSize);
	std::vector<double> u = conjGrad(K, F, u0);

	Coeffs uxCoef{ a1, a2, a3 };
	Coeffs uyCoef{ b1, b2, b3 };

	double normDelta = 0, normU = 0;
	for (int i = 0; i < gSize; ++i)
	{
		double x = grid[i].x;
		double y = grid[i].y;
		std::vector<double> uE = uEx(x, y, uxCoef, uyCoef);
		double ux = u[2 * i];
		double uy = u[2 * i + 1];
		double uN = sqrt(uE[0] * uE[0] + uE[1] * uE[1]);
		if (uN > normU)
			normU = uN;
		double r = sqrt(pow(ux - uE[0], 2) + pow(uy - uE[1], 2));
		if (r > normDelta)
			normDelta = r;
	}
	double delta = normDelta / normU;
	std::cout << delta;

	return 0;
}