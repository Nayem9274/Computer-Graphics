#include <bits/stdc++.h>
#include <cmath>
#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <GL/glut.h>
#endif

using namespace std;

#define pi (2 * acos(0.0))

int drawgrid, drawaxes, clockwise, anticlockwise;
int triangular_face, sphere_segment, cylinder_segment;
float dihedral_angle, cylinder_centre_angle;
float cylinder_h, cylinder_r; // height and radius of cylinder respectively
GLfloat base[][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}};
GLfloat rotate_angle = 0.0f;

float scale = 0.0f;

float sphereColors[][3] = {
    {0, 1, 0}, // Green
    {0, 0, 1}, // Blue
    {0, 1, 0}, // Green
    {0, 0, 1}, // Blue
    {1, 0, 0}, // Red
    {1, 0, 0}  // Red
};

struct Point
{
    double x, y, z;
    Point(double x = 0.0, double y = 0.0, double z = 0.0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Point operator+(Point p)
    {
        Point res;
        res.x = x + p.x;
        res.y = y + p.y;
        res.z = z + p.z;
        return res;
    }
    Point operator*(double k)
    {
        Point res;
        res.x = x * k;
        res.y = y * k;
        res.z = z * k;
        return res;
    }
    Point operator-(Point p)
    {
        Point res;
        res.x = x - p.x;
        res.y = y - p.y;
        res.z = z - p.z;
        return res;
    }
    bool operator==(Point p)
    {
        if ((x == p.x) && (y == p.y) && (z == p.z))
            return true;
        else
        {
            return false;
        }
    }
};

vector<vector<Point>> SpherePoints;

Point center;
Point pos; // camera position
Point l;   // look vector for camera
Point r;   // right vector for camera
Point u;   // up direction for camera
double A = pi / 90;

float Dot_Product(Point u, Point v)
{
    float dot = 0;
    dot += u.x * v.x;
    dot += u.y * v.y;
    dot += u.z * v.z;

    return dot;
}

Point Cross_Product(Point u, Point v)
{
    Point cross;
    cross.x = u.y * v.z - u.z * v.y;
    cross.y = u.z * v.x - u.x * v.z;
    cross.z = u.x * v.y - u.y * v.x;

    return cross;
}

Point Normalize(Point a)
{
    Point result;
    double length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    result.x = a.x / length;
    result.y = a.y / length;
    result.z = a.z / length;

    return result;
}

Point Rotation(Point v, Point reference, int direction)
{
    Point cross = Cross_Product(v, reference); // determining the axis of rotation.
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

float DegreeToRadian(float degree)
{
    return degree * acos(-1.0) / 180.0;
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

void drawTriangle()
{
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < 3; i++)
    {
        glVertex3fv(base[i]);
    }
    glEnd();
}

void drawOctahedron()
{
    // Loop through half of the triangular faces (as they're drawn in pairs)
    for (int i = 0; i < triangular_face / 2; i++)
    {
        // Set alternating colors for the first set of faces
        if (i % 2 == 0)
            glColor3f(1.0f, 0.0f, 1.0f); // Magenta color
        else
            glColor3f(0.0f, 1.0f, 1.0f); // Cyan color

        glPushMatrix();                                      // Save the current transformation state
        glRotatef(i * 90.0, 0, 1, 0);                        // Rotate to position for each triangular face
        glTranslatef(scale / 3.0, scale / 3.0, scale / 3.0); // Translate to the center of the shape along Z-axis
        glScalef(1 - scale, 1 - scale, 1 - scale);           // Scale down the triangle [here 1 is the max scaling]
        drawTriangle();                                      // Draw the first triangle face
        glPopMatrix();                                       // Restore previous transformation state

        // Set alternating colors for the second set of faces
        if (i % 2 == 0)
            glColor3f(0.0f, 1.0f, 1.0f); // Cyan color
        else
            glColor3f(1.0f, 0.0f, 1.0f); // Magenta color

        glPushMatrix();                                      // Save the current transformation state
        glRotatef(i * 90.0, 0, 1, 0);                        // Rotate to position for each triangular face
        glRotatef(180.0, 1, 0, 1);                           // Rotate to draw the second triangle (Orient the second set of triangles symmetric to the first set)
        glTranslatef(scale / 3.0, scale / 3.0, scale / 3.0); // Translate to the center of the shape along Z-axis
        glScalef(1 - scale, 1 - scale, 1 - scale);           // Scale down the triangle
        drawTriangle();                                      // Draw the second triangle face
        glPopMatrix();                                       // Restore previous transformation state
    }
}

void drawCylinder(int segments, float cylinder_r, float cylinder_h)
{
    struct Point coord[segments + 1];
    // Translate and rotate to position the cylinder correctly
    cylinder_h *= (1 - scale);
    cylinder_r *= scale;
    glTranslatef((1 - scale) / sqrt(2), 0, 0);
    glRotatef(-cylinder_centre_angle / 2, 0, 0, 1);
    // Calculate coordinates for each segment of the cylinder
    for (int i = 0; i < segments + 1; i++)
    {
        float theta = i * DegreeToRadian(cylinder_centre_angle) / segments;
        coord[i].x = cylinder_r * cos(theta);
        coord[i].y = cylinder_r * sin(theta);
    }

    glBegin(GL_QUADS);
    for (int i = 0; i < segments; i++)
    {
        //  Drawing each quad for the cylinder
        glVertex3f(coord[i].x, coord[i].y, cylinder_h / 2);
        glVertex3f(coord[i].x, coord[i].y, -cylinder_h / 2);
        glVertex3f(coord[i + 1].x, coord[i + 1].y, -cylinder_h / 2);
        glVertex3f(coord[i + 1].x, coord[i + 1].y, cylinder_h / 2);
    }
    glEnd();
}

void drawAllCylinders()
{
    glColor3f(1.0, 1.0, 0.0);

    // Draw cylinders along the Y-axis rotations
    for (int i = 0; i < 4; i++)
    {
        glPushMatrix();
        glRotatef(45.0, 0, 1, 0);                 // Rotate along Y-axis
        drawCylinder(50, cylinder_r, cylinder_h); // Draw cylinder
        glPopMatrix();

        glRotatef(90.0, 0, 1, 0); // Increment rotation for the next cylinder
    }

    // Draw cylinders along the X-axis rotations
    glPushMatrix();
    glRotatef(90.0, 1, 0, 0); // Rotate along X-axis
    glRotatef(45.0, 0, 1, 0); // Rotate along Y-axis
    for (int i = 0; i < 4; i++)
    {
        glPushMatrix();
        drawCylinder(50, cylinder_r, cylinder_h); // Draw cylinder
        glPopMatrix();

        glRotatef(90.0, 0, 1, 0); // Increment rotation for the next cylinder
    }
    glPopMatrix();

    // Draw cylinders along the Z-axis rotations
    glPushMatrix();
    glRotatef(90.0, 0, 0, 1); // Rotate along Z-axis
    glRotatef(45.0, 0, 1, 0); // Rotate along Y-axis
    for (int i = 0; i < 4; i++)
    {
        glPushMatrix();
        drawCylinder(50, cylinder_r, cylinder_h); // Draw cylinder
        glPopMatrix();

        glRotatef(90.0, 0, 1, 0); // Increment rotation for the next cylinder
    }
    glPopMatrix();
}

// https://www.songho.ca/opengl/gl_sphere.html
// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
vector<vector<Point>> buildUnitPositiveX(int subdivision)
{

    vector<vector<Point>> sphere_points;
    vector<float> vertices;
    float n1[3]; // normal of longitudinal plane rotating along Y-axis
    float n2[3]; // normal of latitudinal plane rotating along Z-axis
    float v[3];  // direction vector intersecting 2 planes, n1 x n2
    float a1;    // longitudinal angle along Y-axis
    float a2;    // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for (unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DegreeToRadian((45.0f - 90.0f * i / (pointsPerRow - 1)));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        vector<Point> points;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for (unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DegreeToRadian((-45.0f + 90.0f * j / (pointsPerRow - 1)));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scaled = 1 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            v[0] *= scaled;
            v[1] *= scaled;
            v[2] *= scaled;

            // add a vertex into array
            vertices.push_back(v[0] / sqrt(3));
            vertices.push_back(v[1] / sqrt(3));
            vertices.push_back(v[2] / sqrt(3));
            points.push_back({v[0] / sqrt(3), v[1] / sqrt(3), v[2] / sqrt(3)});
        }

        sphere_points.push_back(points);
    }

    return sphere_points;
}

void drawSphere(float *colors)
{

    glPushMatrix();
    glColor3fv(colors);
    // Translate and scale the sphere
    glTranslatef(1 - scale, 0.0, 0.0);
    glScalef(scale, scale, scale);
    // Drawing the sphere using QUADS
    glBegin(GL_QUADS);
    for (int k = 0; k < SpherePoints.size() - 1; k++)
    {
        for (int l = 0; l < SpherePoints[l].size() - 1; l++)
        {
            glVertex3f(SpherePoints[k][l].x, SpherePoints[k][l].y, SpherePoints[k][l].z);
            glVertex3f(SpherePoints[k][l + 1].x, SpherePoints[k][l + 1].y, SpherePoints[k][l + 1].z);
            glVertex3f(SpherePoints[k + 1][l + 1].x, SpherePoints[k + 1][l + 1].y, SpherePoints[k + 1][l + 1].z);
            glVertex3f(SpherePoints[k + 1][l].x, SpherePoints[k + 1][l].y, SpherePoints[k + 1][l].z);
        }
    }
    glEnd();
    glPopMatrix();
}

void drawAllSpheres()
{
    glPushMatrix();
    // Draw the upper hemisphere spheres
    for (int i = 0; i <= sphere_segment / 2; i++)
    {
        drawSphere(sphereColors[i]);
        glRotatef(90.0, 0.0, 1.0, 0.0); // Rotate for the next sphere
    }
    // Rotate for the lower hemisphere sphere
    glRotatef(90.0, 0.0, 0.0, 1.0);

    for (int i = sphere_segment / 2 + 1; i < sphere_segment; i++)
    {
        drawSphere(sphereColors[i]);
        glRotatef(180.0, 0.0, 0.0, 1.0);
    }
    glPopMatrix();
}

void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
    case ',':
        // makingthe triangles smaller, spheres larger
        scale += 1.0 / 20.0;
        if (scale >= 1.0)
            scale = 1.0;
        break;
    case '.':
        // making the triangle larger, spheres smaller
        scale -= 1.0 / 20.0;
        if (scale <= 0.0)
            scale = 0.0;
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

    case 'a': // a - rotate the object in the clockwise direction about its own axis
        rotate_angle -= 10.0;
        if (rotate_angle <= -360.0)
            rotate_angle += 360.0;
        break;

    case 'd': // d – rotate the object in the counter-clockwise direction about its own axis.
        rotate_angle += 10.0;
        if (rotate_angle >= -360.0)
            rotate_angle -= 360.0;
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
    glPushMatrix();
    glRotatef(rotate_angle, 0, 1, 0);
    drawOctahedron();
    drawAllCylinders();
    drawAllSpheres();
    glPopMatrix();

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
    drawaxes = 0;

    clockwise = 1;
    anticlockwise = -1;

    pos.x = 4;
    pos.y = 4;
    pos.z = 4;

    l.x = -4;
    l.y = -4;
    l.z = -4;

    r.x = 1;
    r.y = 0;
    r.z = 0;

    u.x = 0;
    u.y = 1;
    u.z = 0;

    triangular_face = 8;
    cylinder_segment = 12;
    sphere_segment = 6;

    dihedral_angle = acos(-1.0 / 3.0) * 180.0 / pi;
    cylinder_centre_angle = 180.0 - dihedral_angle;

    cylinder_h = sqrt(2);
    cylinder_r = 1.0 / sqrt(3);

    SpherePoints = buildUnitPositiveX(3);

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

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshapeListener(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();            // Reset the projection matrix

    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(50, 50);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("Task 3 : Magic Cube ");

    init();

    glEnable(GL_DEPTH_TEST); // enable Depth Testing

    glutDisplayFunc(display); // display callback function
    glutIdleFunc(animate);    // what you want to do in the idle time (when no drawing is occuring)
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop(); // The main loop of OpenGL

    return 0;
}