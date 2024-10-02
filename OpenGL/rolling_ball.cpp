#include <bits/stdc++.h>

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <GL/glut.h>
#endif

using namespace std;

#define pi (2 * acos(0.0))

int drawgrid, drawaxes;
int clockwise = 1;
int anticlockwise = -1;
int wise = 1, go_forward = -1,length_square=200;
bool paused;
double radius, rotation_angle, change_angle, angleX, angleY;
float speedConstant = 0.5;

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
Point coord, direction;
Point center;
Point pos; // camera position
Point l;   // camera ls at
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

    double length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    Point result;
    result.x = a.x / length;
    result.y = a.y / length;
    result.z = a.z / length;

    return result;
}

Point getReflected(Point d, Point n)
{
    Point r;
    // cout<<n.x<<endl<<n.y<<endl<<n.z<<endl;
    double c = (2 * Dot_Product(d, n));
    r = d - n * c;
    return r;
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

void drawGrid()
{
    if (drawgrid == 1)
    {
        float squareSize = 50.0f; // Size of each square

        for (int i = -8; i <= 8; i++)
        {
            for (int j = -8; j <= 8; j++)
            {
                if ((i + j) % 2 == 0)   // Check if the sum of i and j is even for alternating colors
                    glColor3f(1, 1, 1); // Light color for one square
                else
                    glColor3f(0, 0, 0); // Dark color for the other square

                glBegin(GL_QUADS);
                {
                    glVertex3f(i * squareSize, j * squareSize, 0);
                    glVertex3f((i + 1) * squareSize, j * squareSize, 0);
                    glVertex3f((i + 1) * squareSize, (j + 1) * squareSize, 0);
                    glVertex3f(i * squareSize, (j + 1) * squareSize, 0);
                }
                glEnd();
            }
        }
    }
}

const int sectors = 24; // Number of sectors (longitude lines)
const int stacks = 12;  // Number of stacks (latitude lines)

void drawSphere()
{
    float radius = 15.0f;
    float sectorStep = 2.0f * M_PI / sectors;
    float stackStep = M_PI / stacks;
    bool makeTexCoords = true;

    for (int j = 0; j < stacks; j++)
    {
        double latitude1 = (M_PI / stacks) * j - M_PI / 2;
        double latitude2 = (M_PI / stacks) * (j + 1) - M_PI / 2;
        double sinLat1 = sin(latitude1);
        double cosLat1 = cos(latitude1);
        double sinLat2 = sin(latitude2);
        double cosLat2 = cos(latitude2);
        glBegin(GL_QUAD_STRIP);
        for (int i = 0; i <= sectors; i++)
        {
            double longitude = (2 * M_PI / sectors) * i;
            double sinLong = sin(longitude);
            double cosLong = cos(longitude);
            double x1 = cosLong * cosLat1;
            double y1 = sinLong * cosLat1;
            double z1 = sinLat1;
            double x2 = cosLong * cosLat2;
            double y2 = sinLong * cosLat2;
            double z2 = sinLat2;
            // Alternating colors
            if ((i + j) % 2 == 0)
                glColor3f(0.0, 1.0, 0.0); // Green
            else
                glColor3f(1.0, 0.0, 0.0); // Red
            glNormal3d(x2, y2, z2);
            if (makeTexCoords)
                glTexCoord2d(1.0 / sectors * i, 1.0 / stacks * (j + 1));
            glVertex3d(radius * x2, radius * y2, radius * z2);
            glNormal3d(x1, y1, z1);
            if (makeTexCoords)
                glTexCoord2d(1.0 / sectors * i, 1.0 / stacks * j);
            glVertex3d(radius * x1, radius * y1, radius * z1);
        }
        glEnd();
    }
}

void drawArrow()
{
    glPushMatrix();
    {
        // Translate to the surface of the sphere
        // glTranslatef(center.x, center.y+radius, center.z); // Translate to the top surface of the sphere

        // Get the normalized direction vector for the arrow
        // For example, if the ball's direction vector is 'l', use it here
        Point arrowDirection = l; // Change this to the ball's direction

        // Calculate the rotation angles to align the arrow with the direction vector
        angleX = asin(arrowDirection.y);                    // Calculate the angle for X-axis rotation
        angleY = atan2(arrowDirection.x, arrowDirection.z); // Calculate the angle for Y-axis rotation

        // Rotate the arrow to point in the same direction as the ball
        glRotatef(angleX * 180 / pi, 1.0, 0.0, 0.0); // Rotate around the X-axis
        glRotatef(angleY * 180 / pi, 0.0, 1.0, 0.0); // Rotate around the Y-axis

        // glRotatef(rotation_angle, 0.0, 1, 0); // Rotate around the Z-axis by 'rotation_angle' degrees

        // Draw the arrow body (cylinder as shaft)
        glColor3f(0, 255, 255); // Color for the arrow shaft
        GLUquadricObj *cylinder = gluNewQuadric();
        gluCylinder(cylinder, 0.5, 0.5, 25.0, 20, 20); // Parameters: base radius, top radius, height, slices, stacks

        // Draw the arrowhead as a cone
        glTranslatef(0.0, 0.0, 25.0);    // Translate to the tip of the arrow shaft
        glColor3f(0, 0.0, 255);          // Color for the arrowhead
        glutSolidCone(2.0, 6.0, 10, 10); // Draw a cone (base radius, height, slices, stacks)
        glPopMatrix();
    }
    glPopMatrix();
}

void drawSquare(double a)
{
    glLineWidth(10);
    glBegin(GL_LINES);
    {
        glVertex3f(a, a, 0);
        glVertex3f(a, -a, 0);

        glVertex3f(-a, -a, 0);
        glVertex3f(-a, a, 0);

        glVertex3f(-a, -a, 0);
        glVertex3f(a, -a, 0);

        glVertex3f(a, a, 0);
        glVertex3f(-a, a, 0);
    }
    glEnd();
    glLineWidth(1);
}

void keyboardListener(unsigned char key, int x, int y)
{
    double distanceTraveled, angularDisplacement;
    switch (key)
    {

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

    case 'j': // w - move up without changing reference point
        change_angle += 20.0f;
        break;

    case 'l': // s - move down without changing reference point
        change_angle -= 20.0f;
        break;

    case 'w': // w - move up without changing reference point
        pos.y += 0.1;
        l.y -= 0.1;
        break;

    case 's': // s - move down without changing reference point
        pos.y -= 0.1;
        l.y += 0.1; // PRB
        break;

    case ' ':
        paused = !paused;
        break;

    case 'i':
        if (paused)
        {
            go_forward = 1;
            rotation_angle += 45;
            if (rotation_angle >= 360.0)
                rotation_angle -= 360.0;

            coord.x -= speedConstant * radius * pi / 6.0 * cos(change_angle * pi / 180);
            coord.y -= speedConstant * radius * pi / 6.0 * sin(change_angle * pi / 180);
        }

        break;

    case 'k':
        if (paused)
        {
            go_forward = 0;
            rotation_angle -= 45;
            if (rotation_angle < 0.0)
                rotation_angle += 360.0;

            coord.x += speedConstant * radius * pi / 6.0 * cos(change_angle * pi / 180);
            coord.y += speedConstant * radius * pi / 6.0 * sin(change_angle * pi / 180);
        }
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
        pos.x += l.x;
        pos.y += l.y;
        pos.z += l.z;
        break;

    case GLUT_KEY_DOWN: // move backward
        pos.x -= l.x;
        pos.y -= l.y;
        pos.z -= l.z;

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

    drawGrid();
    drawAxes();
    glColor3f(128, 0, 128);
    drawSquare(length_square);

    glPushMatrix();
    {

        // glRotatef(90,1,0,0);

        if (!paused)
        {
            glTranslatef(coord.x, coord.y, coord.z);
            glRotatef(rotation_angle, 1, 0, 0);
            // glRotatef(change_angle, 0, 0, 1);
        }
        else
        {
            glTranslatef(coord.x, coord.y, coord.z);
            // glRotatef(90, 1, 0, 0);
            // glRotatef(change_angle, 0, 1, 0);
            glRotatef(rotation_angle, 1, 0, 0);
        }
        drawSphere();
    }

    glPopMatrix();
    glTranslatef(coord.x, coord.y, coord.z + radius); // Translate to the sphere's position
    glRotatef(change_angle, 0, 0, 1);

    drawArrow();

    // ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate()
{
    if (!paused)
    {
        if (coord.x >= (length_square - radius) && abs(int(coord.x)) > (length_square - radius))
        {
            /*direction = getReflected(direction, Point(-1, 0, 0));
            rotation_angle = -rotation_angle;
            wise = -1;
            cout << "1" << endl;*/
            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(-1, 0, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            rotation_angle = -rotation_angle;
        }
        else if (coord.x <= -(length_square - radius) && abs(int(coord.x)) > (length_square - radius))
        {
            /*direction = getReflected(direction, Point(1, 0, 0));
            rotation_angle = -rotation_angle;
            wise = 1;
            cout << "2" << endl;*/
            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(1, 0, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            rotation_angle = -rotation_angle;
        }
        if (coord.y >= (length_square - radius) && abs(int(coord.y)) > (length_square - radius))
        {
            /* direction = getReflected(direction, Point(0, -1, 0));
              rotation_angle = -rotation_angle;
              wise == -1;
              cout << "3" << endl;*/
            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(0, -1, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            rotation_angle = -rotation_angle;
        }
        else if (coord.y <= -(length_square - radius) && abs(int(coord.y)) > (length_square - radius))
        {
            /*direction = getReflected(direction, Point(0, 1, 0));
            rotation_angle = -rotation_angle;
            wise = 1;
            cout << "4" << endl;*/
            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(0, 1, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            rotation_angle = -rotation_angle;
        }
        // Calculate distance traveled
        // double distanceTraveled = speedConstant * sqrt(direction.x * direction.x + direction.y * direction.y);

        // Calculate angular displacement for rotation
        // double angularDisplacement = distanceTraveled * wise / (2 * pi * radius) * 360; // Convert radians to degrees

        // Update rotation angle
        // rotation_angle += angularDisplacement;
        rotation_angle += 45;

        /* coord.x += speedConstant * direction.x;
         coord.y += speedConstant * direction.y;*/
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        coord.x -= speedConstant * radius * pi / 6.0 * cos(change_angle * pi / 180);
        coord.y -= speedConstant * radius * pi / 6.0 * sin(change_angle * pi / 180);
    }
    else
    {
        // Collision handling logic when paused
        bool collisionHandled = false; // Flag to track if collision was handled

        if (coord.x >= (length_square - radius) && abs(int(coord.x)) > (length_square - radius))
        {

            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(-1, 0, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            // rotation_angle = -rotation_angle;
            collisionHandled = true;
        }
        else if (coord.x <= -(length_square - radius) && abs(int(coord.x)) > (length_square - radius))
        {

            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(1, 0, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            // rotation_angle = -rotation_angle;
            collisionHandled = true;
        }

        if (coord.y >= (length_square - radius) && abs(int(coord.y)) > (length_square - radius))
        {

            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(0, -1, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            //  rotation_angle = -rotation_angle;
            collisionHandled = true;
        }
        else if (coord.y <= -(length_square - radius) && abs(int(coord.y)) > (length_square - radius))
        {

            direction = getReflected(Point(speedConstant * cos(change_angle * pi / 180), speedConstant * sin(change_angle * pi / 180), 0), Point(0, 1, 0));
            change_angle = atan2(direction.y, direction.x) * 180 / pi;
            //  rotation_angle = -rotation_angle;
            collisionHandled = true;
        }
        if (go_forward == 0 && collisionHandled == true)
        {
            coord.x += speedConstant * radius * pi / 6.0 * cos(change_angle * pi / 180);
            coord.y += speedConstant * radius * pi / 6.0 * sin(change_angle * pi / 180);
        }
        else if (go_forward == 1 && collisionHandled == true)
        {

            coord.x -= speedConstant * radius * pi / 6.0 * cos(change_angle * pi / 180);
            coord.y -= speedConstant * radius * pi / 6.0 * sin(change_angle * pi / 180);
        }
    }
    go_forward = -1;

    glutPostRedisplay();
}

void init()
{
    drawgrid = 1;
    drawaxes = 0;
    wise = 1;

    pos.x = 290;
    pos.y = 290;
    pos.z = 90;

    l.x = -1 / sqrt(2.0);
    l.y = -1 / sqrt(2.0);
    l.z = 0;

    r.x = -1 / sqrt(2.0);
    r.y = 1 / sqrt(2.0);
    r.z = 0;

    u.x = 0;
    u.y = 0;
    u.z = 1;

    paused = true;

    // initialize sphere position
    radius = 15.0;
    rotation_angle = 0.0;
    change_angle = 0.0;
    center = Point(21.0, 120.0, radius);

    coord.x = center.x;
    coord.y = center.y;
    coord.z = center.z;
    direction.x = rand() * 0.1 / RAND_MAX;
    direction.y = rand() * 0.1 / RAND_MAX;

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
    gluPerspective(80, 1, 1, 1000.0);
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(50, 50);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("Task 2 : Rolling Ball ");

    init();

    glEnable(GL_DEPTH_TEST); // enable Depth Testing

    glutDisplayFunc(display); // display callback function

    glutIdleFunc(animate); // what you want to do in the idle time (when no drawing is occuring)
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop(); // The main loop of OpenGL

    return 0;
}