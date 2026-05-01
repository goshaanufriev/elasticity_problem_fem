#include <iostream>
#include <vector>
#include <fstream>

struct Node
{
    double x, y;
    bool isBoundary;
};

struct Segment
{
    double a, b;
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
            int k = j * nx + i;
            grid[k].x = xSeg.a + i * hx;
            grid[k].y = ySeg.a + j * hy;
            grid[k].isBoundary = i == 0 || i == nx || j == 0 || j == ny;
        }
    }
    return grid;
}

int Delta(int& i, int& j)
{
    if (i == j)
        return 1;
    return 0;
}

int main()
{
    
    return 0;
}