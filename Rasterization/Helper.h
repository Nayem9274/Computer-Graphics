#include <bits/stdc++.h>
#include "bitmap_image.hpp"
using std::max;
using std::min;
using namespace std;

#define pi (2 * acos(0.0))

static unsigned long int g_seed = 1;

inline int randomInt()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

struct Point
{
    double x, y, z, scale;

    Point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
        this->scale = 1;
    }

    Point(double x, double y, double z) : x(x), y(y), z(z) { scale = 1.0; }

    Point(double x, double y, double z, double scale)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->scale = scale;
    }

    Point(const Point &p)
    {
        this->x = p.x;
        this->y = p.y;
        this->z = p.z;
        this->scale = p.scale;
    }

    void normalize()
    {
        double length = sqrt(x * x + y * y + z * z);
        if (length == 0)
        {
            cout << "Can't divide by zero !!!" << endl;
            return;
        }
        this->x /= length;
        this->y /= length;
        this->z /= length;
    }

    Point operator+(Point p)
    {
        return Point(this->x + p.x, this->y + p.y, this->z + p.z);
    }

    Point operator-(Point p)
    {
        return Point(this->x - p.x, this->y - p.y, this->z - p.z);
    }

    Point operator*(double value)
    {
        return Point(this->x * value, this->y * value, this->z * value);
    }

    Point operator/(double value)
    {
        return Point(this->x / value, this->y / value, this->z / value);
    }

    // DOT PRODUCT
    double dot(Point p)
    {
        return this->x * p.x + this->y * p.y + this->z * p.z;
    }

    // CROSS PRODUCT
    Point operator*(Point p)
    {
        return Point(this->y * p.z - this->z * p.y, this->z * p.x - this->x * p.z, this->x * p.y - this->y * p.x);
    }

    void input(ifstream &fin)
    {
        fin >> x >> y >> z;
        // cout<< x << y << z <<endl;
    }

    void output(ofstream &fout)
    {
        fout << fixed << setprecision(7) << x << " " << y << " " << z << endl;
    }
};

struct Triangle
{
    Point a, b, c;
    rgb_t color;

    Triangle() {}

    Triangle(Point a, Point b, Point c)
    {
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void readTriangle(ifstream &fin)
    {
        a.input(fin);
        b.input(fin);
        c.input(fin);
    }

    void printTriangle(ofstream &fout)
    {
        a.output(fout);
        b.output(fout);
        c.output(fout);
    }

    void scaleTriangle(double sa, double sb, double sc)
    {
        if (sa != 0)
            a = a / sa;
        if (sb != 0)
            b = b / sb;
        if (sc != 0)
            c = c / sc;
    }

    void colorTriangle()
    {
        // static std::mt19937 rng(std::random_device{}());
        // std::uniform_int_distribution<int> dist(0, 255);
        int r = randomInt() % 256;
        int g = randomInt() % 256;
        int b = randomInt() % 256;

        unsigned char red = static_cast<unsigned char>(r);
        unsigned char green = static_cast<unsigned char>(g);
        unsigned char blue = static_cast<unsigned char>(b);
        this->color = make_colour(red, green, blue);
    }

    double getMinX()
    {
        return std::min({a.x, b.x, c.x});
    }

    double getMaxX()
    {
        return std::max({a.x, b.x, c.x});
    }

    double getMinY()
    {
        return std::min({a.y, b.y, c.y});
    }

    double getMaxY()
    {
        return std::max({a.y, b.y, c.y});
    }

    vector<double> getIntersections(double row)
    {
        vector<Point> points = {a, b, c};
        vector<double> intersections;

        for (int i = 0; i < 3; i++)
        {
            int j = (i + 1) % 3;
            if (points[i].y == points[j].y)
                continue;
            if (row >= min(points[i].y, points[j].y) && row <= max(points[i].y, points[j].y))
            {
                intersections.push_back(points[i].x + (row - points[i].y) * (points[i].x - points[j].x) / (points[i].y - points[j].y));
            }
        }
        return intersections;
    }

    // get Z coordinate
    double getZ(double x, double y)
    {
        // find normal of this triangle's plane
        Point A = b - a;
        Point B = c - a;
        Point normal = A * B;
        normal.normalize();
        /* The equation of the plane formed by the triangle is Ax + By + Cz + D = 0, where (A, B, C) are components of the normal vector, and D is a constant.
        Here, D is found using D = -normal.dot(a), where a is a point on the plane (a vertex of the triangle).
        */
        double d = -normal.dot(a);
        // Calculate z-coordinate
        double z = -(normal.x * x + normal.y * y + d) / normal.z;
        return z;
    }
};

struct Matrix
{
    vector<vector<double>> matrix;
    int row; // row=column as square matrix

    Matrix()
    {
        row = 4;
        matrix = vector<vector<double>>(row, vector<double>(row, 0));
    }

    Matrix(int row)
    {
        this->row = row;
        matrix = vector<vector<double>>(row, vector<double>(row, 0));
    }

    Matrix operator*(Matrix m)
    {
        Matrix result(row);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < row; j++)
            {
                for (int k = 0; k < row; k++)
                {
                    result.matrix[i][j] += matrix[i][k] * m.matrix[k][j];
                }
            }
        }
        return result;
    }

    Point operator*(Point p)
    {
        Point result;
        // Assuming row represents the number of rows and columns in the square matrix
        for (int i = 0; i < row; ++i)
        {
            double temp = 0.0;
            for (int j = 0; j < row; ++j)
            {
                temp += matrix[i][j] * ((j == 0) ? p.x : (j == 1) ? p.y
                                                     : (j == 2)   ? p.z
                                                                  : result.scale);
            }
            switch (i)
            {
            case 0:
                result.x = temp;
                break;
            case 1:
                result.y = temp;
                break;
            case 2:
                result.z = temp;
                break;
            case 3:
                result.scale = temp;
                break;
            }
        }
        return result;
    }

    // identity matrix
    void setIdentityMatrix()
    {
        matrix = vector<vector<double>>(row, vector<double>(row, 0));
        for (int i = 0; i < row; i++)
        {
            matrix[i][i] = 1;
        }
    }

    // translation matrix
    void setTranslationMatrix(double tx, double ty, double tz)
    {
        setIdentityMatrix();
        matrix[0][3] = tx;
        matrix[1][3] = ty;
        matrix[2][3] = tz;
    }

    // scaling matrix
    void setScaleMatrix(double sx, double sy, double sz)
    {
        setIdentityMatrix();
        matrix[0][0] = sx;
        matrix[1][1] = sy;
        matrix[2][2] = sz;
    }

    // rodrigues formula
    Point rodriguesFormula(Point x, Point a, double theta)
    {
        double thetaRad = theta * pi / 180.0;
        double cosTheta = cos(thetaRad);
        double sinTheta = sin(thetaRad);
        Point result = x * cosTheta + a * (a.dot(x)) * (1 - cosTheta) + (a * x) * sinTheta;
        return result;
    }

    // rotation matrix
    void setRotationMatrix(double theta, double ax, double ay, double az)
    {
        setIdentityMatrix();
        Point a(ax, ay, az);
        a.normalize();
        Point principalAxes[3];
        principalAxes[0] = Point(1, 0, 0); // i
        principalAxes[1] = Point(0, 1, 0); // j
        principalAxes[2] = Point(0, 0, 1); // k

        for (int i = 0; i < 3; i++)
        {
            Point c = rodriguesFormula(principalAxes[i], a, theta);
            matrix[0][i] = c.x;
            matrix[1][i] = c.y;
            matrix[2][i] = c.z;
        }
    }

    // View Matrix
    void viewMatrix(Point eye, Point look, Point up)
    {
        setIdentityMatrix();
        Point l = look - eye;
        l.normalize();
        Point r = l * up;
        r.normalize();
        Point u = r * l;

        matrix[0][0] = r.x;
        matrix[0][1] = r.y;
        matrix[0][2] = r.z;

        matrix[1][0] = u.x;
        matrix[1][1] = u.y;
        matrix[1][2] = u.z;

        matrix[2][0] = -l.x;
        matrix[2][1] = -l.y;
        matrix[2][2] = -l.z;

        Matrix T;
        T.setTranslationMatrix(-eye.x, -eye.y, -eye.z);

        // V=RT
        *this = *this * T;
    }

    // Projection matrix
    void projectionMatrix(double fovY, double aspectRatio, double near, double far)
    {
        setIdentityMatrix();
        double fovX = fovY * aspectRatio;
        double t = near * tan(fovY * pi / 360.0); // t = near * tan(fovY/2) . why 360(180*2)? dividing by 180 for degree to radian mode.
        double r = near * tan(fovX * pi / 360.0);

        matrix[0][0] = near / r;
        matrix[1][1] = near / t;
        matrix[2][2] = -(far + near) / (far - near);
        matrix[2][3] = -(2.0 * far * near) / (far - near);
        matrix[3][2] = -1.0;
        matrix[3][3] = 0.0;
    }

    void printMatrix(ofstream &fout)
    {
        for (int i = 0; i < row; i++)
        {
            fout << fixed << setprecision(7);
            for (int j = 0; j < row; j++)
            {
                fout << matrix[i][j] << " ";
            }
            fout << endl;
        }
    }
};