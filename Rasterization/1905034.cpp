#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include "Helper.h"
using namespace std;

#define pi (2 * acos(0.0))

Point eye, look, up;
double fovY, aspectRatio, near, far;
Matrix M;
stack<Matrix> S;

int total_triangles;

void modelingTransformation()
{
    ifstream fin("../Offline_2/IOs/5/scene.txt");
    ofstream fout("stage1.txt");

    eye.input(fin);
    look.input(fin);
    up.input(fin);

    fin >> fovY >> aspectRatio >> near >> far;
    // cout<<fovY << aspectRatio <<near <<far<<endl;

    M.setIdentityMatrix();
    S.push(M);

    string command;

    while (true)
    {
        fin >> command;
        // cout << command << endl;
        if (command == "triangle")
        {
            Triangle tri;
            tri.readTriangle(fin);
            tri.a = S.top() * tri.a;
            tri.b = S.top() * tri.b;
            tri.c = S.top() * tri.c;

            tri.printTriangle(fout);

            fout << endl;

            total_triangles++;
        }
        else if (command == "translate")
        {
            double tx, ty, tz;
            fin >> tx >> ty >> tz;
            Matrix T;
            T.setTranslationMatrix(tx, ty, tz);
            M = S.top() * T;
            S.pop();
            S.push(M);
        }
        else if (command == "scale")
        {
            double sx, sy, sz;
            fin >> sx >> sy >> sz;
            Matrix T;
            T.setScaleMatrix(sx, sy, sz);
            M = S.top() * T;
            S.pop();
            S.push(M);
        }
        else if (command == "rotate")
        {
            double angle, ax, ay, az;
            fin >> angle >> ax >> ay >> az;
            Matrix T;
            T.setRotationMatrix(angle, ax, ay, az);
            M = S.top() * T;
            S.pop();
            S.push(M);
        }
        else if (command == "push")
        {
            S.push(S.top());
        }
        else if (command == "pop")
        {
            if (S.empty())
            {
                cout << "Stack is empty, nothing to pop" << endl;
                break;
            }
            S.pop();
        }
        else if (command == "end")
        {
            break;
        }
        else
        {
            cout << "Invalid Command! " << command << endl;
            break;
        }
    }
    fin.close();
    fout.close();
}

void viewTransformation()
{
    ifstream fin("stage1.txt");
    ofstream fout("stage2.txt");
    Matrix V;
    V.viewMatrix(eye, look, up);
    for (int i = 0; i < total_triangles; i++)
    {
        Triangle tri;
        tri.readTriangle(fin);

        tri.a = V * tri.a;
        tri.b = V * tri.b;
        tri.c = V * tri.c;

        tri.printTriangle(fout);
        fout << endl;
    }

    fin.close();
    fout.close();
}

void projectionTransformation()
{
    ifstream fin("stage2.txt");
    ofstream fout("stage3.txt");
    Matrix P;
    // cout << fovY << aspectRatio << near << far << endl;
    P.projectionMatrix(fovY, aspectRatio, near, far);
    for (int i = 0; i < total_triangles; i++)
    {
        Triangle tri;
        tri.readTriangle(fin);

        tri.a = P * tri.a;
        tri.b = P * tri.b;
        tri.c = P * tri.c;

        tri.scaleTriangle(tri.a.scale, tri.b.scale, tri.c.scale);
        tri.printTriangle(fout);

        fout << endl;
    }

    fin.close();
    fout.close();
}

void clippingScanConversion()
{
    // 1. Read data
    ifstream config("../Offline_2/IOs/5/config.txt");
    int screenWidth, screenHeight;
    config >> screenWidth >> screenHeight;
    config.close();

    ifstream fin("stage3.txt");
    Triangle triangles[total_triangles];
    for (int i = 0; i < total_triangles; i++)
    {
        triangles[i].readTriangle(fin);
        triangles[i].colorTriangle();
    }
    fin.close();

    // 2. Initialize z buffer and frame buffer
    double leftLimit = -1.0, rightLimit = 1.0; // X-axis
    double bottomLimit = -1.0, topLimit = 1.0; // Y- axis
    double zMax = 1.0, zMin = -1.0;

    double dx = (rightLimit - leftLimit) / screenWidth;
    double dy = (topLimit - bottomLimit) / screenHeight;

    double topY = topLimit - (dy / 2.0);
    double bottomY = bottomLimit + (dy / 2.0);
    double leftX = leftLimit + (dx / 2.0);
    double rightX = rightLimit - (dx / 2.0);

    vector<vector<double>> zBuffer(screenHeight, vector<double>(screenWidth, zMax));

    bitmap_image image(screenWidth, screenHeight);
    for (int i = 0; i < screenWidth; i++)
    {
        for (int j = 0; j < screenHeight; j++)
        {
            image.set_pixel(i, j, 0, 0, 0); // initialize background color with black
        }
    }

    // 3. Apply z-buffer algorithm
    for (int i = 0; i < total_triangles; i++)
    {
        double min_x = triangles[i].getMinX(), max_x = triangles[i].getMaxX();
        double min_y = triangles[i].getMinY(), max_y = triangles[i].getMaxY();

        double left_scanline = max(leftX, min_x);
        double right_scanline = min(rightX, max_x);
        double bottom_scanline = max(bottomY, min_y); // bottom scan line
        double top_scanline = min(topY, max_y);       // top scan line

        for (double row = top_scanline; row >= bottom_scanline; row -= dy)
        {
            vector<double> intersections = triangles[i].getIntersections(row);

            if (intersections[0] > intersections[1])
            {
                swap(intersections[0], intersections[1]);
            }

            // Clipping against left and right edges
            intersections[0] = max(intersections[0], left_scanline);
            intersections[1] = min(intersections[1], right_scanline);

            for (double col = intersections[0]; col <= intersections[1]; col += dx)
            {
                double z = triangles[i].getZ(col, row); // Calculate z values
                int colIndex = (col - leftLimit) / dx;
                int rowIndex = (topLimit - row) / dy;
                //Compare with z-buffer and z_front_limit(zMin) and update if required 
                if (z < zBuffer[rowIndex][colIndex] && z > zMin)
                {
                    zBuffer[rowIndex][colIndex] = z;
                    image.set_pixel(colIndex, rowIndex, triangles[i].color);
                }
            }
        }
    }

    // save the z-buffer values( are less than z_max) in z-buffer.txt
    ofstream fout("z_buffer.txt");

    for (int row = 0; row < screenHeight; row++)
    {
        for (int col = 0; col < screenWidth; col++)
        {
            if (zBuffer[row][col] < zMax)
            {
                fout << fixed << setprecision(6) << zBuffer[row][col] << "\t";
            }
        }
        fout << endl;
    }

    fout.close();

    image.save_image("out.bmp");

    // free memory
    zBuffer.clear();
    // cout << "Done !! " << endl;
}

int main()
{
    modelingTransformation();
    viewTransformation();
    projectionTransformation();
    clippingScanConversion();
}
