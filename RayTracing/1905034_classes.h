#include <bits/stdc++.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2 * acos(0.0))

class Object;
static unsigned long int g_seed = 1;

inline int randomInt()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

/* Taken from Helper.h of Offline_2*/
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

    double length() { return sqrt(x * x + y * y + z * z); } // Extra added in this offline

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

    double operator*(Point b) { return x * b.x + y * b.y + z * b.z; }

    // DOT PRODUCT
    double dot(Point p)
    {
        return this->x * p.x + this->y * p.y + this->z * p.z;
    }

    // CROSS PRODUCT
    Point operator&(Point p)
    {
        return Point(this->y * p.z - this->z * p.y, this->z * p.x - this->x * p.z, this->x * p.y - this->y * p.x);
    }

    Point operator-()
    {
        return Point(-x, -y, -z); // Extra added
    }

    void input(ifstream &fin)
    {
        fin >> x >> y >> z;
        // cout<< x << y << z <<endl;
    }

    // Extra added
    friend ifstream &operator>>(ifstream &fin, Point &p)
    {
        fin >> p.x >> p.y >> p.z;
        return fin;
    }

    void output(ofstream &fout)
    {
        fout << fixed << setprecision(7) << x << " " << y << " " << z << endl;
    }
};

// Ray
struct Ray
{
    Point start, dir;

    Ray() {}
    Ray(Point start, Point dir)
    {
        this->start = start;
        dir.normalize();
        this->dir = dir;
    }
};

// Color in RGB
struct Color
{
    double r, g, b;
    Color()
    {
        this->r = 0;
        this->g = 0;
        this->b = 0;
    }
    Color(double r, double g, double b) : r(r), g(g), b(b) {}

    Color operator+(Color c)
    {
        return Color(this->r + c.r, this->g + c.g, this->b + c.b);
    }

    Color operator-(Color c)
    {
        return Color(this->r - c.r, this->g - c.g, this->b - c.b);
    }

    Color operator*(Color c)
    {
        return Color(this->r * c.r, this->g * c.g, this->b * c.b);
    }

    Color operator+(double value)
    {
        return Color(this->r + value, this->g + value, this->b + value);
    }

    Color operator-(double value)
    {
        return Color(this->r - value, this->g - value, this->b - value);
    }

    Color operator*(double value)
    {
        return Color(this->r * value, this->g * value, this->b * value);
    }

    void operator=(const Color &c)
    {
        this->r = c.r;
        this->g = c.g;
        this->b = c.b;
    }
};

void drawCheckerBoard(double floorWidth, double tileWidth, Point reference_point)
{
    int x = 0;
    int y = 0;
    // Calculate the number of rows and columns based on floorWidth and tileWidth
    int rows = round(floorWidth / tileWidth);
    int cols = round(floorWidth / tileWidth);

    glBegin(GL_QUADS);
    for (int i = x - rows; i < x + rows; i++)
    {
        for (int j = y - cols; j < y + cols; j++)
        {
            // Determining the color of the current tile based on the sum of row and column indices
            if ((i + j) % 2 == 0)
            {
                glColor3f(1.0, 1.0, 1.0); // white
            }
            else
            {
                glColor3f(0.0, 0.0, 0.0); // black
            }
            // // Drawing the current tile as a quadrilateral (4 vertices)
            glVertex3f(reference_point.x + i * tileWidth, reference_point.y + j * tileWidth, 0);
            glVertex3f(reference_point.x + i * tileWidth + tileWidth, reference_point.y + j * tileWidth, 0);
            glVertex3f(reference_point.x + i * tileWidth + tileWidth, reference_point.y + j * tileWidth + tileWidth, 0);
            glVertex3f(reference_point.x + i * tileWidth, reference_point.y + j * tileWidth + tileWidth, 0);
        }
    }
    glEnd();
}

// Determinant calculation for Barycentric triangle calculation
double determinant(double a[3][3])
{
    double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    return det;
}

// Draw Sphere
void drawSphere(double radius, int stacks, int slices, Color color)
{
    glColor3d(color.r, color.g, color.b);
    struct Point points[stacks + 1][slices + 1]; // an array to store the points of the sphere
    for (int j = 0; j <= stacks; j++)
    {
        double phi = -M_PI / 2.0 + j * M_PI / stacks; // Calculate the polar angle phi
        double r = radius * cos(phi);                 // Calculate the radius of the sphere at the current phi
        double h = radius * sin(phi);                 // Calculate the height of the sphere at the current phi

        // Calculate the points on the circle at the current phi
        for (int i = 0; i < slices + 1; i++)
        {
            double theta = i * 2.0 * M_PI / slices;             // Calculate the azimuthal angle theta
            points[j][i] = {r * cos(theta), r * sin(theta), h}; // Calculate the x, y, and z coordinates of the point
        }
    }

    glBegin(GL_QUADS);
    for (int j = 0; j < stacks; j++)
    {
        for (int i = 0; i < slices; i++)
        {
            glVertex3f(points[j][i].x, points[j][i].y, points[j][i].z);
            glVertex3f(points[j][i + 1].x, points[j][i + 1].y, points[j][i + 1].z);
            glVertex3f(points[j + 1][i + 1].x, points[j + 1][i + 1].y, points[j + 1][i + 1].z);
            glVertex3f(points[j + 1][i].x, points[j + 1][i].y, points[j + 1][i].z);
        }
    }
    glEnd();
}

// Draw Cone : Taken from https://www.freemancw.com/2012/06/opengl-cone-function/ & Offline-1 cylinder function
void drawCone(double height, double radius, int segments, Color color)
{
    double tempx = radius, tempy = 0;
    double currx, curry;
    glBegin(GL_TRIANGLES);
    for (int i = 1; i <= segments; i++)
    {
        double theta = i * 2.0 * M_PI / segments;
        currx = radius * cos(theta);
        curry = radius * sin(theta);

        GLfloat c = (2 + cos(theta)) / 3;
        glColor3f(color.r, color.g, color.b);
        glVertex3f(0, 0, height / 2);          // Apex
        glVertex3f(currx, curry, -height / 2); // Base
        glVertex3f(tempx, tempy, -height / 2); // Base
        // Update for next iteration
        tempx = currx;
        tempy = curry;
    }
    glEnd();
}

// Pointlight
struct PointLight
{
    Point position;
    Color color;

    void draw()
    {

        glPushMatrix();
        glTranslated(position.x, position.y, position.z);
        drawSphere(1, 20, 20, {color.r, color.g, color.b});
        glPopMatrix();
    }

    void input(ifstream &fin)
    {
        fin >> position.x >> position.y >> position.z;
        fin >> color.r >> color.g >> color.b;
        if (fin.fail())
        {
            cout << "Error in input taking for PointLight !" << endl;
        }
        // cout<< x << y << z <<endl;
    }

    friend ifstream &operator>>(ifstream &fin, PointLight &l)
    {
        fin >> l.position.x >> l.position.y >> l.position.z;
        fin >> l.color.r >> l.color.g >> l.color.b;
        return fin;
    }
};

// Spotlight
struct SpotLight
{
    PointLight point_light;
    Point light_direction;
    double cutoffAngle;

    void draw()
    {
        Color color = point_light.color;
        Point position = point_light.position;
        glPushMatrix();
        glTranslated(position.x, position.y, position.z);
        // glRotated(-acos(light_direction.x) * 90.0 / M_PI, light_direction.y, light_direction.z, 0.0);
        Point up(0, 0, 1);
        Point axis = up & (light_direction);
        double angle = acos(up.dot(light_direction));
        glRotated((angle + 180) * 180 / (acos(-1.0)), axis.x, axis.y, axis.z);
        drawCone(10, 4, 50, {color.r, color.g, color.b});
        glPopMatrix();
    }

    void input(ifstream &fin)
    {
        point_light.position.input(fin);
        fin >> point_light.color.r >> point_light.color.g >> point_light.color.b;
        light_direction.input(fin);
        fin >> cutoffAngle;
        if (fin.fail())
        {
            cout << "Error in input taking for SpotLight !" << endl;
        }
        // cout<< x << y << z <<endl;
    }

    friend ifstream &operator>>(ifstream &fin, SpotLight &l)
    {
        fin >> l.point_light.position;
        fin >> l.point_light.color.r >> l.point_light.color.g >> l.point_light.color.b;
        fin >> l.light_direction;
        fin >> l.cutoffAngle;
        return fin;
    }
};

// declaration
extern vector<Object *> objects;
extern vector<PointLight *> pointLights;
extern vector<SpotLight *> spotLights;
extern int recursion_level;

// Object
class Object
{
public:
    Point reference_point; // should have x,y,z
    double height, width, length;
    Color color;
    vector<double> coefficients; // ambient, diffuse, specular, reflection coefficients
    int shine;                   // exponent term of specular component

    Object()
    {
        coefficients = vector<double>(4, 0);
        shine = 0;
    }

    virtual void draw() = 0;
    virtual Ray getNormal(Point point, Ray incidentRay) = 0;
    virtual Color getColorAt(Point point)
    {
        return Color(this->color.r, this->color.g, this->color.b);
    }

    virtual double intersectObject(Ray ray, Color &color, int level) = 0; // Implementation different for different objects

    void setColor(Color color)
    {
        this->color = color;
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoefficients(vector<double> coefficients)
    {
        this->coefficients = coefficients;
    }

    void handlePointLight(Point intersectionPoint, const Color &intersectionPointColor, Color &color, Ray ray)
    {
        for (int i = 0; i < pointLights.size(); i++)
        {
            Point light_pos = pointLights[i]->position;
            Point lightDirection = intersectionPoint - light_pos;
            lightDirection.normalize();

            // cast rayl from pl.light_pos to intersectionPoint
            Ray lightRay = Ray(light_pos, lightDirection);

            // calculate normal at intersectionPoint
            Ray normal = getNormal(intersectionPoint, lightRay);

            // Find reflected ray : formula taken from slide W09 pg23
            Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir * 2 * (lightRay.dir * normal.dir));

            // Check if incident ray is not obstructed by any object
            double t2 = (intersectionPoint - light_pos).length();
            if (t2 < 1e-7)
                continue;

            bool obscured = false;
            // checking if the light ray is obstructed by any objects in the scene.
            for (Object *object : objects)
            {
                double t3 = object->intersectObject(lightRay, color, 0);
                /* If the light ray intersects an object and the intersection point is between the light source and the surface point,
                   the light is considered to be obstructed.*/
                if (t3 > 0 && t3 + 1e-7 < t2)
                {
                    obscured = true;
                    break;
                }
            }

            if (!obscured)
            {
                // Lambert value: max {(L.N),0} , -ve sign for correcting direction [see slide W09 : pg21]
                double lambert_value = max(0.0, -lightRay.dir.dot(normal.dir));

                // Here ray is the viewing direction(eye) : {(R.V)^k ,0}]
                double phong_value = max(0.0, -ray.dir.dot(reflection.dir));

                // Update diffuse and specular components. Ambient already updated
                // I = Iaka   + ∑(i=1 to n) Ipi {kd (L.N) + ks (R.V)^k } , k-shining factor
                color.r += pointLights[i]->color.r * coefficients[1] * lambert_value * intersectionPointColor.r;
                color.r += pointLights[i]->color.r * coefficients[2] * pow(phong_value, shine) * intersectionPointColor.r;

                color.g += pointLights[i]->color.g * coefficients[1] * lambert_value * intersectionPointColor.g;
                color.g += pointLights[i]->color.g * coefficients[2] * pow(phong_value, shine) * intersectionPointColor.g;

                color.b += pointLights[i]->color.b * coefficients[1] * lambert_value * intersectionPointColor.b;
                color.b += pointLights[i]->color.b * coefficients[2] * pow(phong_value, shine) * intersectionPointColor.b;
            }
        }
    }

    void handleSpotlight(Point intersectionPoint, const Color &intersectionPointColor, Color &color, Ray ray)
    {
        for (int i = 0; i < spotLights.size(); i++)
        {
            Point light_pos = spotLights[i]->point_light.position;
            Point lightDirection = intersectionPoint - light_pos;
            lightDirection.normalize();

            // Do the same calculation for each spot light unless
            // the ray cast from light_pos to intersectionPoint
            // exceeds cutoff-angle for the light source
            double dotProduct = lightDirection.dot(spotLights[i]->light_direction); // Beta = Dot(A,B)/|A||B| , A-spotlight dir, B-lightDir
            // angle is the angle between spotlight direction and ray casted from spotlight to intersection point. See slide W09:pg-9
            double Betaangle = acos(dotProduct / (lightDirection.length() * spotLights[i]->light_direction.length())) * (180.0 / pi);

            if (fabs(Betaangle) < spotLights[i]->cutoffAngle)
            {
                Ray lightRay = Ray(light_pos, lightDirection);
                Ray normal = getNormal(intersectionPoint, lightRay);

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir * 2 * (lightRay.dir * normal.dir));

                double t2 = (intersectionPoint - light_pos).length();
                if (t2 < 1e-7)
                    continue;

                bool obscured = false;

                for (Object *object : objects)
                {
                    double t3 = object->intersectObject(lightRay, color, 0);
                    if (t3 > 0 && t3 + 1e-7 < t2)
                    {
                        obscured = true;
                        break;
                    }
                }

                if (!obscured)
                {
                    double lambert_value = max(0.0, -lightRay.dir.dot(normal.dir));
                    double phong_value = max(0.0, -ray.dir.dot(reflection.dir));

                    color.r += spotLights[i]->point_light.color.r * coefficients[1] * lambert_value * intersectionPointColor.r;
                    color.r += spotLights[i]->point_light.color.r * coefficients[2] * pow(phong_value, shine) * intersectionPointColor.r;

                    color.g += spotLights[i]->point_light.color.g * coefficients[1] * lambert_value * intersectionPointColor.g;
                    color.g += spotLights[i]->point_light.color.g * coefficients[2] * pow(phong_value, shine) * intersectionPointColor.g;

                    color.b += spotLights[i]->point_light.color.b * coefficients[1] * lambert_value * intersectionPointColor.b;
                    color.b += spotLights[i]->point_light.color.b * coefficients[2] * pow(phong_value, shine) * intersectionPointColor.b;
                }
            }
        }
    }
    void handleRecursiveReflection(Point intersectionPoint, Ray ray, Color &color, int level)
    {
        if (level < recursion_level)
        {
            Ray normal = getNormal(intersectionPoint, ray);
            // construct reflected ray from intersection point
            Ray reflectedRay = Ray(intersectionPoint, ray.dir - normal.dir * 2 * (ray.dir * normal.dir));
            // actually slightly forward from the point (by moving the start a little bit towards the reflection direction) to avoid self intersection
            reflectedRay.start = reflectedRay.start + reflectedRay.dir * 1e-7;

            int nearest = -1; // index of nearest object 
            double t = -1, tMin = 1e7;
            /* find tmin from the nearest intersecting object, using intersect() method, as done in the capture() method
            if found, call intersect(rreflected, colorreflected, level+1) */
            for (int k = 0; k < objects.size(); k++)
            {
                t = objects[k]->intersect(reflectedRay, color, 0);
                if (t > 0 && t < tMin)
                    tMin = t, nearest = k;
            }
            // if nearest object is found, then calculate the color of the object at the point of intersection
            if (nearest != -1)
            {
                Color colorReflected(0, 0, 0);
                double t = objects[nearest]->intersect(reflectedRay, colorReflected, level + 1);
                color.r += colorReflected.r * coefficients[3];
                color.g += colorReflected.g * coefficients[3];
                color.b += colorReflected.b * coefficients[3];
            }
        }
    }

    // ray- ray casted from eye. also it is V(viewing direction) for specular reflection calculation
    virtual double intersect(Ray ray, Color &color, int level)
    {
        // code for finding intersecting tmin
        double t = intersectObject(ray, color, level);
        if (t < 0)
            return -1;
        //if level is 0, return tmin    
        if (level == 0)
            return t;

        Point intersectionPoint = ray.start + ray.dir * t;
        Color intersectionPointColor = getColorAt(intersectionPoint);

        color.r = intersectionPointColor.r * coefficients[0]; // 0- ambience
        color.g = intersectionPointColor.g * coefficients[0];
        color.b = intersectionPointColor.b * coefficients[0];

        handlePointLight(intersectionPoint, intersectionPointColor, color, ray);
        handleSpotlight(intersectionPoint, intersectionPointColor, color, ray);
        handleRecursiveReflection(intersectionPoint, ray, color, level);

        return t;
    }
    // Destructor
    virtual ~Object()
    {
        coefficients.clear();
    }
};

// Triangle
struct Triangle : public Object
{
    Point a, b, c, n; // 3 points and normal

    Triangle() {}

    Triangle(Point a, Point b, Point c)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->n = (b - a) & (c - a);
        this->n.normalize();
    }

    friend ifstream &operator>>(ifstream &fin, Triangle &t)
    {
        fin >> t.a >> t.b >> t.c;                   // x,y,z
        fin >> t.color.r >> t.color.g >> t.color.b; // color
        for (int i = 0; i < 4; i++)
            fin >> t.coefficients[i]; // ambient, diffuse, specular, recursive reflection coefficient
        fin >> t.shine;               // shininess
        t.n = (t.b - t.a) & (t.c - t.a);
        t.n.normalize();
        return fin;
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

    virtual Ray getNormal(Point point, Ray incidentRay)
    {
        if (incidentRay.dir.dot(n) < 0)
        {
            return Ray(point, -n); // the incident ray is coming from the inside of the surface
        }
        else
        {
            return Ray(point, n); // the incident ray is coming from the outside of the surface.
        }
    }

    virtual void draw()
    {
        glColor3f(color.r, color.g, color.b);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    // Intersection with Barycentric Triangle
    //  See slide W10: pg 49,50,51 for detailed formula and condition
    virtual double intersectObject(Ray ray, Color &color, int level)
    {

        double detBeta[3][3] = {
            {a.x - ray.start.x, a.x - c.x, ray.dir.x},
            {a.y - ray.start.y, a.y - c.y, ray.dir.y},
            {a.z - ray.start.z, a.z - c.z, ray.dir.z}};
        double detGamma[3][3] = {
            {a.x - b.x, a.x - ray.start.x, ray.dir.x},
            {a.y - b.y, a.y - ray.start.y, ray.dir.y},
            {a.z - b.z, a.z - ray.start.z, ray.dir.z}};
        double detT[3][3] = {
            {a.x - b.x, a.x - c.x, a.x - ray.start.x},
            {a.y - b.y, a.y - c.y, a.y - ray.start.y},
            {a.z - b.z, a.z - c.z, a.z - ray.start.z}};
        double matA[3][3]{
            {a.x - b.x, a.x - c.x, ray.dir.x},
            {a.y - b.y, a.y - c.y, ray.dir.y},
            {a.z - b.z, a.z - c.z, ray.dir.z}};

        double detA = determinant(matA);
        double Beta = determinant(detBeta) / detA;
        double Gamma = determinant(detGamma) / detA;
        double t = determinant(detT) / detA;

        if (Beta + Gamma < 1 && Beta > 0 && Gamma > 0 && t > 0)
        {
            return t;
        }
        else
        {
            return -1;
        }
    }
};

// Sphere
struct Sphere : public Object
{
    Point center;
    double radius;

    Sphere() {}
    Sphere(Point center, double radius)
    {
        this->center = center;
        this->radius = radius;
    }

    friend ifstream &operator>>(ifstream &fin, Sphere &s)
    {
        fin >> s.center >> s.radius;                // center and radius
        fin >> s.color.r >> s.color.g >> s.color.b; // color
        for (int i = 0; i < 4; i++)
            fin >> s.coefficients[i]; // ambient, diffuse, specular, recursive reflection coefficient
        fin >> s.shine;               // shininess
        return fin;
    }

    /*The normal vector at any point on the surface of a sphere is simply the vector from the center of the sphere to that point on the surface
    In the getNormal function, start is the point of intersection (point), and dir is the vector from the center of the sphere to the point of
    intersection (point - center). So, point - center gives the vector from the center of the sphere to the point of intersection.*/
    virtual Ray getNormal(Point point, Ray incidentRay)
    {
        return Ray(point, point - center);
    }

    virtual void draw()
    {
        glPushMatrix();
        glTranslated(center.x, center.y, center.z);
        drawSphere(radius, 50, 50, color);
        glPopMatrix();
    }

    // Slide W10:Pg 30-31(Ray-Sphere Intersection)
    virtual double intersectObject(Ray ray, Color &color, int level)
    {
        ray.start = ray.start - center;

        double a = 1; // Normalized
        double b = 2 * (ray.dir.dot(ray.start));
        double c = (ray.start.dot(ray.start)) - (radius * radius);

        double discriminant = b * b - 4.0 * a * c;

        if (discriminant < 0)
        {
            return -1;
        }
        // Add later: if abs(a) close to 1e-7 i.e 0  t=-c/b  [Not very important]

        double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
        double t2 = (-b + sqrt(discriminant)) / (2.0 * a);

        return t1 > 0.0 ? t1 : t2;
    }
};

struct Floor : public Object
{
    double tileWidth, floorWidth;
    int tiles;
    Floor() {}
    Floor(double floorWidth, double tileWidth)
    {
        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;
        tiles = round(floorWidth / tileWidth);
        reference_point = Point(-floorWidth / 2, -floorWidth / 2, 0);
        length = tileWidth;
    }

    virtual Ray getNormal(Point point, Ray incidentRay)
    {
        if (incidentRay.dir.z > 0)
            return Ray(point, Point(0, 0, 1)); // If the incidentRay is coming from above the floor, normal vector points upwards from the floor's surface.
        else
            return Ray(point, Point(0, 0, -1)); // If the incidentRay is coming from below the floor, normal vector points downwards from the floor's surface.
    }

    virtual void draw()
    {
        glPushMatrix();
        drawCheckerBoard(floorWidth, tileWidth, reference_point);
        glPopMatrix();
    }

    virtual Color getColorAt(Point point)
    {
        int X = floor((point.x - reference_point.x) / tileWidth); // gives the tile coordinates in terms of the number of tiles from the reference point.
        int Y = floor((point.y - reference_point.y) / tileWidth);

        if ((abs(X) + abs(Y)) % 2 == 0)
        {
            return Color(1, 1, 1);
        }
        else
        {
            return Color(0, 0, 0);
        }
    }

    // See slide W10: pg 24
    virtual double intersectObject(Ray ray, Color &color, int level)
    {
        Point normal = Point(0, 0, 1); // As per offline_information page-10 (Indicates a plane in xy plane)
        double D = 0;                  // Passes through the start
        // checking if the dot product is close to zero within a tolerance of 0.001
        if (round(normal.dot(ray.dir) * 1000) == 0)
            return -1;
        double t = -(D + normal.dot(ray.start)) / (normal.dot(ray.dir));
        // The calculated point p is checked against the bounds of the checkerboard. If it falls outside the bounds, the function returns -1 to indicate no intersection.
        Point p = ray.start + ray.dir * t;
        if (p.x <= reference_point.x || p.x >= abs(reference_point.x) && p.y <= reference_point.y && p.y >= abs(reference_point.y))
        {
            return -1;
        }
        return t;
    }
};

struct GeneralQuadraticSurface : public Object
{
    double A, B, C, D, E, F, G, H, I, J;

    GeneralQuadraticSurface() {}

    friend ifstream &operator>>(ifstream &fin, GeneralQuadraticSurface &g)
    {
        fin >> g.A >> g.B >> g.C >> g.D >> g.E >> g.F >> g.G >> g.H >> g.I >> g.J;
        fin >> g.reference_point >> g.length >> g.width >> g.height;

        fin >> g.color.r >> g.color.g >> g.color.b;
        for (int i = 0; i < 4; i++)
            fin >> g.coefficients[i];
        fin >> g.shine;
        return fin;
    }

    virtual Ray getNormal(Point point, Ray incidentRay)
    {
        // ∂F/∂x = 2Ax + Dy + Ez + G, ∂F/∂y = 2By + Dx + Fz + H, ∂F/∂z = 2Cz + Ex + Fy + I --> Normal(𝜕F/𝜕x, 𝜕F/𝜕y, 𝜕F/𝜕z)
        Point normal(2 * A * point.x + D * point.y + E * point.z + G,
                     2 * B * point.y + D * point.x + F * point.z + H,
                     2 * C * point.z + E * point.x + F * point.y + I);

        return Ray(point, normal);
    }

    virtual void draw()
    {
        ; // you do not have to draw the general quadric surfaces
    }

    bool clipping(Point point)
    {
        // Check if the point lies within the reference cube along the x-axis
        if (length > 0 && (point.x < reference_point.x || point.x > reference_point.x + length))
            return false;

        // Check if the point lies within the reference cube along the y-axis
        if (width > 0 && (point.y < reference_point.y || point.y > reference_point.y + width))
            return false;

        // Check if the point lies within the reference cube along the z-axis
        if (height > 0 && (point.z < reference_point.z || point.z > reference_point.z + height))
            return false;

        // If the point lies within the reference cube along all dimensions, return true
        return true;
    }

    // See slide for detailed calculation Ray Casting & Tracing_OhioState.pdf [pg-14,15] & specs [pg -10]
    virtual double intersectObject(Ray ray, Color &color, int level)
    {
        double Xr = ray.start.x;
        double Yr = ray.start.y;
        double Zr = ray.start.z;

        double Xd = ray.dir.x;
        double Yd = ray.dir.y;
        double Zd = ray.dir.z;

        double a = A * Xd * Xd + B * Yd * Yd + C * Zd * Zd + D * Xd * Yd + E * Xd * Zd + F * Yd * Zd;
        double b = 2 * A * Xr * Xd + 2 * B * Yr * Yd + 2 * C * Zr * Zd + D * (Xr * Yd + Xd * Yr) + E * (Xr * Zd + Xd * Zr) + F * (Yr * Zd + Yd * Zr) + G * Xd + H * Yd + I * Zd;
        double c = A * Xr * Xr + B * Yr * Yr + C * Zr * Zr + D * Xr * Yr + E * Xr * Zr + F * Yr * Zr + G * Xr + H * Yr + I * Zr + J;

        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0)
            return -1;
        // if the quadratic equation has a linear solution. i.e a close to zero then b*t+c=0 ---> t=-c/b
        if (fabs(a) < 1e-7)
        {
            return -c / b;
        }
        double t1 = (-b - sqrt(discriminant)) / (2 * a);
        double t2 = (-b + sqrt(discriminant)) / (2 * a);

        if (t1 < 0 && t2 < 0)
            return -1;

        if (t2 < t1)
            swap(t1, t2); // Swap t1 and t2 so that t1 is the smaller positive solution

        // Check if t1 is positive and within bounds, if so, return t1
        if (t1 > 0)
        {
            Point intersectionPoint = ray.start + ray.dir * t1;
            if (clipping(intersectionPoint))
            {
                return t1;
            }
        }
        //  Check if t2 is positive and within bounds, if so, return t2
        if (t2 > 0)
        {
            Point intersectionPoint = ray.start + ray.dir * t2;
            if (clipping(intersectionPoint))
            {
                return t2;
            }
        }

        return -1;
    }
};
