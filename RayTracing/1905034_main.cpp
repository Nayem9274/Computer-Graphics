#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include "1905034_classes.h"
#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <GL/glut.h>
#endif

using namespace std;

#define pi (2 * acos(0.0))

int drawgrid, drawaxes, clockwise, anticlockwise;
double A = pi / 90;
Point pos; // camera position
Point l;   // look vector for camera
Point r;   // right vector for camera
Point u;   // up direction for camera

int recursion_level, objectCount, pointlightCount, spotlightCount;
int imageHeight, imageWidth; // same
bitmap_image image;
int imageCount = 0;
int windowWidth = 640;
int windowHeight = 640;
double viewAngle = 80.0;

vector<Object *> objects;
vector<PointLight *> pointLights;
vector<SpotLight *> spotLights;

// Taken from Offline-1
Point Rotation(Point v, Point reference, int direction)
{
    Point cross = v & reference; // determining the axis of rotation.
    Point output;

    double B = direction * A;
    /*The cos(A) term usually represents the component of the original vector along the axis of rotation that remains unchanged,
    while the sin(B) term, multiplied by the result of the cross product between the original vector v and the axis of rotation (cross),
    represents the contribution to the new vector due to rotation.*/

    output.x = v.x * cos(A) + cross.x * sin(B);
    output.y = v.y * cos(A) + cross.y * sin(B);
    output.z = v.z * cos(A) + cross.z * sin(B);

    return output;
}

void capture()
{
    cout << "Bitmap Image Capturing........." << endl;

    bitmap_image image(imageWidth, imageHeight); // initialize bitmap image
    for (int column = 0; column < imageWidth; column++)
    {
        for (int row = 0; row < imageHeight; row++)
        {
            image.set_pixel(column, row, 0, 0, 0); // set background color
        }
    }
    double planeDistance = windowHeight / (2.0 * tan(viewAngle / 2.0 * pi / 180.0));
    Point topLeft = pos + l * planeDistance - r * (windowWidth / 2.0) + u * (windowHeight / 2.0);

    double du = ((double)windowWidth / imageWidth);
    double dv = ((double)windowHeight / imageHeight);

    // Choose middle of the grid cell
    topLeft = topLeft + r * (du / 2.0) - u * (dv / 2.0);
    int nearest = -1;
    double t, tMin;
    for (int i = 0; i < imageWidth; i++)
    {
        for (int j = 0; j < imageHeight; j++)
        {
            // calculate curPixel using topleft,r,u,i,j,du,dv
            Point curPixel = topLeft + (r * du * i) - (u * dv * j);

            // cast ray from eye i.e camera(pos) to (curPixel-eye) direction
            Ray ray(pos, curPixel - pos);
            Color color;
            // save the nearest object
            tMin = -1;
            nearest = -1; // index of nearest object 
            for (int k = 0; k < objects.size(); k++)
            {
                t = objects[k]->intersect(ray, color, 0);
                if (t > 0 && (nearest == -1 || t < tMin))
                    tMin = t, nearest = k;
            }

            // if nearest object is found, then calculate the color of the object at the point of intersection
            if (nearest != -1)
            {

                color = Color(0, 0, 0);
                double t = objects[nearest]->intersect(ray, color, 1);
                // convert the floating-point RGB values to integer RGB values 
                int red = round(color.r * 255.0);
                int green = round(color.g * 255.0);
                int blue = round(color.b * 255.0);
                image.set_pixel(i, j, red, green, blue);
            }
        }
    }
    /*  save image // The 1st image you capture after running the program should be named Output_11.bmp, the 2nd image you capture should be
          named Output_12.bmp and so on.*/
    imageCount++;
    image.save_image("Output_1" + to_string(imageCount) + ".bmp");
    cout << "Bitmap Image Captured!!!" << endl;
}

void drawAxes()
{
    if (drawaxes == 1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f(100, 0, 0);
            glVertex3f(-100, 0, 0);

            glVertex3f(0, -100, 0);
            glVertex3f(0, 100, 0);

            glVertex3f(0, 0, 100);
            glVertex3f(0, 0, -100);
        }
        glEnd();
    }
}

void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '0':
        capture();
        break;

    case '1': // Look left
        l = Rotation(l, u, anticlockwise);
        r = Rotation(r, u, anticlockwise);
        break;

    case '2': // Look right
        l = Rotation(l, u, clockwise);
        r = Rotation(r, u, clockwise);
        break;

    case '3': // Look up
        l = Rotation(l, r, anticlockwise);
        u = Rotation(u, r, anticlockwise);
        break;

    case '4': // Look down
        l = Rotation(l, r, clockwise);
        u = Rotation(u, r, clockwise);
        break;

    case '5': // tilt clockwise
        r = Rotation(r, l, anticlockwise);
        u = Rotation(u, l, anticlockwise);
        break;

    case '6': // tilt counter clockwise
        r = Rotation(r, l, clockwise);
        u = Rotation(u, l, clockwise);
        break;

    case 'w': // w - move up without changing reference point
        pos.y += 0.1;
        l.y -= 0.1; // PRB

        break;

    case 's': // s - move down without changing reference point
        pos.y -= 0.1;
        l.y += 0.1; // PRB

        break;

    default:
        break;
    }
}

void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP: // move forward
        pos.x += 0.6 * l.x;
        pos.y += 0.6 * l.y;
        pos.z += 0.6 * l.z;
        break;

    case GLUT_KEY_DOWN: // move backward
        pos.x -= 0.6 * l.x;
        pos.y -= 0.6 * l.y;
        pos.z -= 0.6 * l.z;
        break;

    case GLUT_KEY_RIGHT: //  move right
        pos.x += r.x;
        pos.y += r.y;
        pos.z += r.z;
        break;

    case GLUT_KEY_LEFT: // move left
        pos.x -= r.x;
        pos.y -= r.y;
        pos.z -= r.z;
        break;

    case GLUT_KEY_PAGE_UP: // page up, move up
        pos.x += u.x;
        pos.y += u.y;
        pos.z += u.z;
        break;

    case GLUT_KEY_PAGE_DOWN: // page down, move down
        pos.x -= u.x;
        pos.y -= u.y;
        pos.z -= u.z;
        break;
    default:
        break;
    }
}

void mouseListener(int button, int state, int x, int y)
{
}

void loadData()
{
    ifstream fin("scene.txt");
    fin >> recursion_level >> imageHeight;
    imageWidth = imageHeight;
    fin >> objectCount;
    for (int i = 0; i < objectCount; i++)
    {
        string object;
        fin >> object;
        Object *obj;

        if (object == "sphere")
        {
            obj = new Sphere();
            fin >> *((Sphere *)obj);
        }
        else if (object == "triangle")
        {
            obj = new Triangle();
            fin >> *((Triangle *)obj);
        }
        else if (object == "general")
        {
            obj = new GeneralQuadraticSurface();
            fin >> *((GeneralQuadraticSurface *)obj);
        }
        else
        {
            cout << object << " is not a valid object!!" << endl;
        }
        objects.push_back(obj);
    }

    fin >> pointlightCount;

    for (int i = 0; i < pointlightCount; i++)
    {
        PointLight *light = new PointLight();
        fin >> *light;
        pointLights.push_back(light);
    }

    fin >> spotlightCount;

    for (int i = 0; i < spotlightCount; i++)
    {
        SpotLight *spotlight = new SpotLight();
        fin >> *spotlight;
        spotLights.push_back(spotlight);
    }

    Object *floor;
    floor = new Floor(1000, 20);
    floor->setColor(Color(0.85, 0.85, 0.85));
    vector<double> coefficients = {0.4, 0.2, 0.2, 0.2};
    floor->setCoefficients(coefficients);
    objects.push_back(floor);
}

void display()
{
    // clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0); // color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-u camera here
    ********************/
    // load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    // initialize the matrix
    glLoadIdentity();
    // control viewing (or camera)

    // camera control
    gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);

    // again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    drawAxes();

    for (int i = 0; i < objects.size(); i++)
    {
        objects[i]->draw();
    }

    for (int i = 0; i < pointLights.size(); i++)
    {
        pointLights[i]->draw();
    }

    for (int i = 0; i < spotLights.size(); i++)
    {
        spotLights[i]->draw();
    }

    // ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate()
{
    glutPostRedisplay();
}

void init()
{
    drawgrid = 1;
    loadData();
    drawaxes = 0;

    clockwise = 1;
    anticlockwise = -1;

    pos.x = 100,
    pos.y = 100,
    pos.z = 10;

    u.x = 0,
    u.y = 0,
    u.z = 1;

    r.x = -1/ sqrt(2),
    r.y = 1/ sqrt(2),
    r.z = 0;

    l.x = -1/ sqrt(2),
    l.y = -1/ sqrt(2),
    l.z = 0;

    // clear the screen
    glClearColor(0, 0, 0, 0);

    /************************
    / set-u projection here
    ************************/
    // load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    // initialize the matrix
    glLoadIdentity();

    // give PERSPECTIVE parameters
    gluPerspective(80, 1, 1, 1000.0);
}



int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("Offline-3 (1905034) : Ray Tracing ");

    init();

    glEnable(GL_DEPTH_TEST); // enable Depth Testing

    glutDisplayFunc(display); // display callback function
    glutIdleFunc(animate);    // what you want to do in the idle time (when no drawing is occuring)
    //glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop(); // The main loop of OpenGL

    // CLEAR MEMORY
    objects.clear();
    pointLights.clear();
    spotLights.clear();

    return 0;
}