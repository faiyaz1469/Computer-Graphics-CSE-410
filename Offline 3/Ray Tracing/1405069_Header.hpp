#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<limits>

using namespace std;

#define PI 2*acos(0.0)
#define INF numeric_limits<double>::infinity()

double _rdn = (PI * 3)/180; // 3 degree,in rad format

struct point
{
	double x,y,z;

	point(){}

	point(double a,double b,double c){
	    x = a;
	    y = b;
	    z = c;
	}

    point normalize(){
        point tmp;
        double magnitude = sqrt(x*x + y*y + z*z);
        tmp.x = x / magnitude;
        tmp.y = y / magnitude;
        tmp.z = z / magnitude;
        return tmp;
    }

    point operator+(const point & temp)const{
        point tmp;
        tmp.x = x + temp.x;
        tmp.y = y + temp.y;
        tmp.z = z + temp.z;
        return tmp;
    }

	point operator-(const point & temp)const{
        point tmp;
        tmp.x = x - temp.x;
        tmp.y = y - temp.y;
        tmp.z = z - temp.z;
        return tmp;
    }

    point operator*(const double & temp)const{
        point tmp;
        tmp.x = x * temp;
        tmp.y = y * temp;
        tmp.z = z * temp;
        return tmp;
    }

    point operator/(const double & temp)const{
        point tmp;
        tmp.x = x / temp;
        tmp.y = y / temp;
        tmp.z = z / temp;
        return tmp;
    }

	point nonvector_product(const point & temp)const{
        point tmp;
        tmp.x = x * temp.x;
        tmp.y = y * temp.y;
        tmp.z = z * temp.z;
        return tmp;
    }

	double dot_product(const point & temp)const{
		double dotProd = x * temp.x + y * temp.y + z * temp.z;
        return dotProd;
    }

    point cross_product(const point & temp)const{
        point tmp;
        tmp.x = y * temp.z - z * temp.y;
        tmp.y = z * temp.x - x * temp.z;
        tmp.z = x * temp.y - y * temp.x;
        return tmp;
    }

	point angle_product(const point & temp)const{
        point tmp;
        tmp.x = x * cos(_rdn) + temp.x * sin(_rdn);
        tmp.y = y * cos(_rdn) + temp.y * sin(_rdn);
        tmp.z = z * cos(_rdn) + temp.z * sin(_rdn);
        return tmp;
    }

	point angle_product_2(const point & temp)const{
        point tmp;
        tmp.x = x * cos(-_rdn) + temp.x * sin(-_rdn);
        tmp.y = y * cos(-_rdn) + temp.y * sin(-_rdn);
        tmp.z = z * cos(-_rdn) + temp.z * sin(-_rdn);
        return tmp;
    }

};


class Vector {
    double x;
    double y;
    double z;

public:
    Vector() {
        x = 0.0;
		y = 0.0;
		z = 0.0;
    }

    Vector(double x1, double y1, double z1) {
        x = x1;
        y = y1;
        z = z1;
    }

	 void setX(double x1) {
        x = x1;
    }

    void setY(double y1) {
        y = y1;
    }

    void setZ(double z1) {
        z = z1;
    }

    double getX() const {
        return x;
    }

    double getY() const {
        return y;
    }

    double getZ() const {
        return z;
    }

    void normalize()
    {
        double r = sqrt((x * x) + (y * y) + (z * z));
        x = x / r;
        y = y / r;
        z = z / r;
    }

    Vector operator+(const Vector vec) {
		double a = x + vec.x;
		double b = y + vec.y;
		double c = z + vec.z;
		Vector tmp(a, b, c);
		return tmp;
	}

    Vector operator-(const Vector vec){
		double a = x - vec.x;
		double b = y - vec.y;
		double c = z - vec.z;
		Vector tmp(a, b, c);
		return tmp;
	}

	double calcDist(Vector vec) {
		double dist = sqrt((x - vec.x)*(x - vec.x) + (y - vec.y)*(y - vec.y) + (z - vec.z)*(z - vec.z));
		return dist;
	}

	Vector scalarProduct(const double scalar){
		double a = x * scalar;
		double b = y * scalar;
		double c = z * scalar;
		Vector tmp(a, b, c);
		return tmp;
	}

	double dotProduct(const Vector vec){
		double a = x * vec.x;
		double b = y * vec.y;
		double c = z * vec.z;
		double dotProd = a + b + c;
		return dotProd;
	}

	Vector crossProduct(const Vector vec){
		double a = y * vec.z - z * vec.y;
		double b = z * vec.x - x * vec.z;
		double c = x * vec.y - y * vec.x;
		Vector tmp(a, b, c);
		return tmp;
	}

    ~Vector() {}
};


class Ray {
    Vector rStart;
    Vector rDirection;

public:
    Ray() {}

    Ray(Vector rOrg, Vector rDir) {
        rStart = rOrg;
        rDirection = rDir;
        rDirection.normalize();
    }

    Vector getRayStart() const {
        return rStart;
    }

    Vector getRayDir() const {
        return rDirection;
    }

    ~Ray() {}
};


struct Color {

	double red;
    double green;
    double blue;

    Color() {
        red = 0.0;
		green = 0.0;
		blue = 0.0;
    }

    Color(double r, double g, double b) {
        red = r;
        green = g;
        blue = b;
    }

    ~Color() {}
};


class PointLight {
    Vector position;
    Color color;

    double radius;
    int sgmnts;
    int stcks;

public:
    PointLight() {
        radius = 0.0;
        sgmnts = 0;
		stcks = 0;
    }

    PointLight(Vector pos, Color clr, double rad, int seg, int stk) {
        position = pos;
        color = clr;

        radius = rad;
        sgmnts = seg;
        stcks = stk;
    }

	 PointLight(Vector pos, Color clr) {
        position = pos;
        color = clr;
    }

    Vector getPosition() const {
        return position;
    }

    Color getColor() const {
        return color;
    }

    void draw();

    ~PointLight() {}
};


void PointLight::draw() {
    Vector pts[stcks+1][sgmnts+1];
    double height, tempRadius;

    // generating points: sgmnts = segments in plane; stcks = segments in hemisphere
	for(int i=0; i<=stcks; i++) {
		double angle = ((double)i/(double)stcks)*(PI/2);
		height = radius * sin(angle);
		tempRadius = radius * cos(angle);

		for(int j=0; j<=sgmnts; j++) {
			double ang = ((double)j/(double)sgmnts)*2*PI;
			double a = tempRadius*cos(ang);
			double b = tempRadius*sin(ang);
			double c = height;
            pts[i][j] = Vector(a, b, c);
		}
	}

	// drawing quads using generated points
	glColor3f(color.red, color.green, color.blue);

	for(int i=0; i<stcks; i++) {
		for(int j=0; j<sgmnts; j++) {
			glBegin(GL_QUADS);
			{
			    Vector newPos1 = position + pts[i][j];
				Vector newPos2 = position + pts[i][j+1];
				Vector newPos3 = position + pts[i+1][j+1];
				Vector newPos4 = position + pts[i+1][j];

				Vector newPos5 = position - pts[i][j];
				Vector newPos6 = position - pts[i][j+1];
				Vector newPos7 = position - pts[i+1][j+1];
				Vector newPos8 = position - pts[i+1][j];

				/* upper hemisphere */
				glVertex3f(newPos1.getX(), newPos1.getY(), newPos1.getZ());
				glVertex3f(newPos2.getX(), newPos2.getY(), newPos2.getZ());
				glVertex3f(newPos3.getX(), newPos3.getY(), newPos3.getZ());
				glVertex3f(newPos4.getX(), newPos4.getY(), newPos4.getZ());

                /* lower hemisphere */
                glVertex3f(newPos1.getX(), newPos1.getY(), newPos5.getZ());
				glVertex3f(newPos2.getX(), newPos2.getY(), newPos6.getZ());
				glVertex3f(newPos3.getX(), newPos3.getY(), newPos7.getZ());
				glVertex3f(newPos4.getX(), newPos4.getY(), newPos8.getZ());

			}
			glEnd();
		}
	}
}

class SpotLight {
	PointLight point_light;
	Vector direction;
	double cutOff_angle;

    double radius;
    int sgmnts;
    int stcks;

public:
    SpotLight() {
        radius = 0.0;
        sgmnts = 0;
		stcks = 0;
    }

    SpotLight(PointLight pLight, Vector dir, double cutAngle, double rad, int seg, int stk) {

		point_light = pLight;
		direction = dir;
		cutOff_angle = cutAngle;

        radius = rad;
        sgmnts = seg;
        stcks = stk;
    }

	Vector getDirection() const {
		return direction;
	}

	double getCutoffAngle() const {
		return cutOff_angle;
	}

	PointLight getPointLight() const {
		return point_light;
	}

    void draw();

    ~SpotLight() {}
};


void SpotLight::draw() {
    Vector pts[stcks+1][sgmnts+1];
    double height, tempRadius;

    // generating points: sgmnts = segments in plane; stcks = segments in hemisphere
	for(int i=0; i<=stcks; i++) {
		double angle = ((double)i/(double)stcks)*(PI/2);
		height = radius * sin(angle);
		tempRadius = radius * cos(angle);

		for(int j=0; j<=sgmnts; j++) {
			double ang = ((double)j/(double)sgmnts)*2*PI;
			double a = tempRadius*cos(ang);
			double b = tempRadius*sin(ang);
			double c = height;
            pts[i][j] = Vector(a, b, c);
		}
	}

	// drawing quads using generated points
	Color color = point_light.getColor();
	glColor3f(color.red, color.green, color.blue);

	Vector position = point_light.getPosition();

	for(int i=0; i<stcks; i++) {
		for(int j=0; j<sgmnts; j++) {
			glBegin(GL_QUADS);
			{
			    Vector newPos1 = position + pts[i][j];
				Vector newPos2 = position + pts[i][j+1];
				Vector newPos3 = position + pts[i+1][j+1];
				Vector newPos4 = position + pts[i+1][j];

				Vector newPos5 = position - pts[i][j];
				Vector newPos6 = position - pts[i][j+1];
				Vector newPos7 = position - pts[i+1][j+1];
				Vector newPos8 = position - pts[i+1][j];

				/* upper hemisphere */
				glVertex3f(newPos1.getX(), newPos1.getY(), newPos1.getZ());
				glVertex3f(newPos2.getX(), newPos2.getY(), newPos2.getZ());
				glVertex3f(newPos3.getX(), newPos3.getY(), newPos3.getZ());
				glVertex3f(newPos4.getX(), newPos4.getY(), newPos4.getZ());

                /* lower hemisphere */
                glVertex3f(newPos1.getX(), newPos1.getY(), newPos5.getZ());
				glVertex3f(newPos2.getX(), newPos2.getY(), newPos6.getZ());
				glVertex3f(newPos3.getX(), newPos3.getY(), newPos7.getZ());
				glVertex3f(newPos4.getX(), newPos4.getY(), newPos8.getZ());
			}
			glEnd();
		}
	}
}


struct ReflectionCoefficient {
    double ambientRC;
    double diffuseRC;
    double specularRC;
    double recursiveRC;

    ReflectionCoefficient() {
        ambientRC = 0.0;
		diffuseRC = 0.0;
		specularRC = 0.0;
		recursiveRC = 0.0;
    }

    ReflectionCoefficient(double amb, double diff, double spec, double rec) {
        ambientRC = amb;
        diffuseRC = diff;
        specularRC = spec;
        recursiveRC = rec;
    }

    ~ReflectionCoefficient() {}
};


class Object {
    Color color;
    ReflectionCoefficient reflectionCoefficient;
    int shininess;

public:
    Object() {
        shininess = 0;
    }

    Color getColor() const {
        return color;
    }

    void setColor(Color clr) {
        color = clr;
    }

    ReflectionCoefficient getReflectionCoefficient() const {
        return reflectionCoefficient;
    }

    void setReflectionCoefficient(ReflectionCoefficient rc) {
        reflectionCoefficient = rc;
    }

    int getShininess() const {
        return shininess;
    }

    void setShininess(int shine) {
        shininess = shine;
    }

    void calcAmbient(Color& color, Color intersectClr) {
		color.red = intersectClr.red * getReflectionCoefficient().ambientRC;
		color.green = intersectClr.green * getReflectionCoefficient().ambientRC;
		color.blue = intersectClr.blue * getReflectionCoefficient().ambientRC;
	}

    void calcReflection(Ray ray, Color& color, Vector intersectPnt, Color intersectClr, Vector normal, PointLight light, Ray incidentRay){

		double lambertVal = (incidentRay.getRayDir().scalarProduct(-1.0)).dotProduct(normal);
		double tempVal = ((incidentRay.getRayDir().dotProduct(normal))*2.0);
		Ray reflectedRay(intersectPnt, incidentRay.getRayDir() - normal.scalarProduct(tempVal));
		double phongVal = (ray.getRayDir().scalarProduct(-1.0)).dotProduct(reflectedRay.getRayDir());
		double fixed = getReflectionCoefficient().diffuseRC * max(lambertVal, 0.0);

		double tempRed = light.getColor().red * intersectClr.red;
		double tempGreen = light.getColor().green * intersectClr.green;
		double tempBlue = light.getColor().blue * intersectClr.blue;

		color.red = color.red + tempRed * fixed;
		color.green = color.green + tempGreen * fixed;
		color.blue = color.blue + tempBlue * fixed;

		double fixedShine = getReflectionCoefficient().specularRC * pow(max(phongVal, 0.0), getShininess());

		color.red = color.red + tempRed * fixedShine;
		color.green = color.green + tempGreen * fixedShine;
		color.blue = color.blue + tempBlue * fixedShine;
	}

    virtual void draw() = 0;
    virtual double intersect(Ray, Color&, int) = 0;

    virtual ~Object() {}
};


Vector position;
int recursionLevel = 0;

vector<Object*> objects;
vector<PointLight> point_lights;
vector<SpotLight> spot_lights;


class Sphere: public Object {
    Vector reference_point;
    double radius;

    int sgmnts;
    int stcks;

public:
    Sphere() {
        radius = 0.0;
        sgmnts = 0;
		stcks = 0;
    }

    Sphere(Vector center, double rad, int seg, int stks) {
        reference_point = center;
        radius = rad;

        sgmnts = seg;
        stcks = stks;
    }

    void draw();
    double intersect(Ray, Color&, int);

    ~Sphere() {}
};

void Sphere::draw() {
    Vector pts[stcks+1][sgmnts+1];
    double height, tempRadius;

	// generating points: sgmnts = segments in plane; stcks = segments in hemisphere
	for(int i=0; i<=stcks; i++) {
		double angle = ((double)i/(double)stcks)*(PI/2);
		height = radius * sin(angle);
		tempRadius = radius * cos(angle);

		for(int j=0; j<=sgmnts; j++) {
			double ang = ((double)j/(double)sgmnts)*2*PI;
			double a = tempRadius*cos(ang);
			double b = tempRadius*sin(ang);
			double c = height;
            pts[i][j] = Vector(a, b, c);
		}
	}

	// drawing quads using generated points
	glColor3f(getColor().red, getColor().green, getColor().blue);

	for(int i=0; i<stcks; i++) {
		for(int j=0; j<sgmnts; j++) {
			glBegin(GL_QUADS);
			{
			    Vector newPos1 = reference_point + pts[i][j];
				Vector newPos2 = reference_point + pts[i][j+1];
				Vector newPos3 = reference_point + pts[i+1][j+1];
				Vector newPos4 = reference_point + pts[i+1][j];

				Vector newPos5 = reference_point - pts[i][j];
				Vector newPos6 = reference_point - pts[i][j+1];
				Vector newPos7 = reference_point - pts[i+1][j+1];
				Vector newPos8 = reference_point - pts[i+1][j];

				/* upper hemisphere */
				glVertex3f(newPos1.getX(), newPos1.getY(), newPos1.getZ());
				glVertex3f(newPos2.getX(), newPos2.getY(), newPos2.getZ());
				glVertex3f(newPos3.getX(), newPos3.getY(), newPos3.getZ());
				glVertex3f(newPos4.getX(), newPos4.getY(), newPos4.getZ());

                /* lower hemisphere */
                glVertex3f(newPos1.getX(), newPos1.getY(), newPos5.getZ());
				glVertex3f(newPos2.getX(), newPos2.getY(), newPos6.getZ());
				glVertex3f(newPos3.getX(), newPos3.getY(), newPos7.getZ());
				glVertex3f(newPos4.getX(), newPos4.getY(), newPos8.getZ());

			}
			glEnd();
		}
	}
}


double Sphere::intersect(Ray ray, Color& color, int level) {

	// finding intersecting tMin
    double a, b, c, discriminant, tMax, tMin, tTemp, tTemp2;

    a = ray.getRayDir().dotProduct(ray.getRayDir());
    b = ((ray.getRayStart().dotProduct(ray.getRayDir())) - (ray.getRayDir().dotProduct(reference_point))) * 2.0;
    c = (ray.getRayStart().dotProduct(ray.getRayStart())) + (reference_point.dotProduct(reference_point)) - (ray.getRayStart().dotProduct(reference_point)) * 2.0 - radius*radius;

    discriminant = b*b-4.0*a*c;
	tTemp = -b/(2.0*a);
	tTemp2 = sqrt(discriminant)/(2.0*a);

    if(discriminant < 0.0) {
        tMin = INF;
    }
	else if(discriminant > 0.0) {
        tMax = tTemp + tTemp2;
        tMin = tTemp - tTemp2;

		if(tMin > 0.0)
			tMin = tMin;
		else
			tMin = tMax;
    }
	else {
        tMin = tTemp;
    }

    if(level == 0) {
        return tMin;
    }

    // illuminating with Phong Lighting Model
    Vector intersectPnt = ray.getRayStart()+ ray.getRayDir().scalarProduct(tMin);
    Color intersectClr = getColor();

    // determining unit normal vector at intersection point on object's surface
    Vector normal = intersectPnt - reference_point;
    normal.normalize();

    if(position.calcDist(reference_point) > radius){
	   normal = normal;
    }
    else{
	   normal = normal.scalarProduct(-1.0);
    }

    // calculating ambient light component of reflected ray
    calcAmbient(color, intersectClr);

    // calculating diffuse & specular reflection components of reflected ray
    for(unsigned int i=0; i<point_lights.size(); i++) {
        Ray incidentRay(point_lights[i].getPosition(), intersectPnt - point_lights[i].getPosition());

        // checking if intersection point is in shadow
        double tmp;
		double tMinimum = INF;

        for(unsigned int j=0; j<objects.size(); j++) {
            Color dummyColor;  // black
            tmp = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(tmp>0.0 && tmp<tMinimum) {
                tMinimum = tmp;
            }
        }

        Vector shadowIntersectionPoint = incidentRay.getRayStart() + incidentRay.getRayDir().scalarProduct(tMinimum);
        double epsilon = 0.0000001;  // for tuning light effect

        if(intersectPnt.calcDist(incidentRay.getRayStart())-epsilon > shadowIntersectionPoint.calcDist(incidentRay.getRayStart())) {
            // intersection point is, indeed, in shadow
            continue;
        }

        calcReflection(ray, color, intersectPnt, intersectClr, normal, point_lights[i], incidentRay);
    }

	if(level >= recursionLevel) {
        return tMin;
    }

    return tMin;
}



class Triangle: public Object {
    Vector a;
    Vector b;
    Vector c;

public:
    Triangle() {}

    Triangle(Vector a1, Vector b1, Vector c1) {
        a = a1;
        b = b1;
        c = c1;
    }

    void draw();
    double intersect(Ray, Color&, int);

    ~Triangle() {}
};

void Triangle::draw() {
    glColor3f(getColor().red, getColor().green, getColor().blue);

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(a.getX(), a.getY(), a.getZ());
        glVertex3f(b.getX(), b.getY(), b.getZ());
        glVertex3f(c.getX(), c.getY(), c.getZ());
    }
    glEnd();
}


double Triangle::intersect(Ray ray, Color& color, int level) {
	// finding intersecting tMin
    double detBase, detBeta, detGamma, detT, tMin;
	double tempBase1, tempBase2, tempBase3, tempBase4, tempBase5, tempBase6, tempBase7, tempBase8, tempBase9;

	tempBase1 = ((a.getY()-c.getY())*ray.getRayDir().getZ()-(a.getZ()-c.getZ())*ray.getRayDir().getY());
	tempBase2 = ((a.getZ()-b.getZ())*ray.getRayDir().getY()-(a.getY()-b.getY())*ray.getRayDir().getZ());
	tempBase3 = ((a.getY()-b.getY())*(a.getZ()-c.getZ())-(a.getZ()-b.getZ())*(a.getY()-c.getY()));

    detBase = (a.getX() - b.getX()) * tempBase1;
    detBase = detBase + (a.getX() - c.getX()) * tempBase2;
    detBase = detBase + ray.getRayDir().getX() * tempBase3;

	tempBase4 = ((a.getZ()-ray.getRayStart().getZ())*ray.getRayDir().getY()-(a.getY()-ray.getRayStart().getY())*ray.getRayDir().getZ());
	tempBase5 = ((a.getY()-ray.getRayStart().getY())*(a.getZ()-c.getZ())-(a.getZ()-ray.getRayStart().getZ())*(a.getY()-c.getY()));

    detBeta = (a.getX() - ray.getRayStart().getX()) * tempBase1;
    detBeta = detBeta + (a.getX() - c.getX()) * tempBase4;
    detBeta = detBeta + ray.getRayDir().getX() * tempBase5;

	tempBase6 = ((a.getY()-ray.getRayStart().getY())*ray.getRayDir().getZ()-(a.getZ()-ray.getRayStart().getZ())*ray.getRayDir().getY());
	tempBase7 = ((a.getY()-b.getY())*(a.getZ()-ray.getRayStart().getZ())-(a.getZ()-b.getZ())*(a.getY()-ray.getRayStart().getY()));

    detGamma = (a.getX() - b.getX()) * tempBase6;
    detGamma = detGamma + (a.getX() - ray.getRayStart().getX()) * tempBase2;
    detGamma = detGamma + ray.getRayDir().getX() * tempBase7;

	tempBase8 = ((a.getY()-c.getY())*(a.getZ()-ray.getRayStart().getZ())-(a.getZ()-c.getZ())*(a.getY()-ray.getRayStart().getY()));
	tempBase9 = ((a.getZ()-b.getZ())*(a.getY()-ray.getRayStart().getY())-(a.getY()-b.getY())*(a.getZ()-ray.getRayStart().getZ()));

    detT = (a.getX() - b.getX()) * tempBase8;
    detT = detT + (a.getX() - c.getX()) * tempBase9;
    detT = detT + (a.getX() - ray.getRayStart().getX()) * tempBase3;

    if(detBase == 0.0) {
        // ray will not intersect the triangle plane
        tMin = INF;
    }
	else {
        // ray will intersect the triangle plane
        if(detBeta/detBase > 0.0 && detGamma/detBase > 0.0 && detBeta/detBase+detGamma/detBase < 1.0) {
            // intersection point lies within the boundary of the triangle
            tMin = detT/detBase;
        } else {
            // intersection point does not lie within the boundary of the triangle
            tMin = INF;
        }
    }

    if(level == 0) {
        return tMin;
    }

    // illuminating with Phong Lighting Model
    Vector intersectPnt = ray.getRayStart() + ray.getRayDir().scalarProduct(tMin);
    Color intersectClr = getColor();

    // determining unit normal vector on appropriate side of triangle
    Vector tmp1 = b-a;
    Vector tmp2 = c-a;
    Vector normal = tmp1.crossProduct(tmp2);
    normal.normalize();

	if((ray.getRayDir().scalarProduct(-1.0)).dotProduct(normal) > 0.0){
	   normal = normal;
    }
    else{
	   normal = normal.scalarProduct(-1.0);
    }

    // calculating ambient light component of reflected ray
    calcAmbient(color, intersectClr);

    // calculating diffuse & specular reflection components of reflected ray
    for(unsigned int i=0; i<point_lights.size(); i++) {
        Ray incidentRay(point_lights[i].getPosition(), intersectPnt-point_lights[i].getPosition());

        // checking if intersection point is in shadow
        double tmp, tMinimum = INF;

        for(unsigned int j=0; j<objects.size(); j++) {
            Color dummyColor;  // black
            tmp = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(tmp>0.0 && tmp<tMinimum) {
                tMinimum = tmp;
            }
        }

        Vector shadowIntersectionPoint = incidentRay.getRayStart() + incidentRay.getRayDir().scalarProduct(tMinimum);
        double epsilon = 0.0000001;  // for tuning light effect

        if(intersectPnt.calcDist(incidentRay.getRayStart())-epsilon > shadowIntersectionPoint.calcDist(incidentRay.getRayStart())) {
            // intersection point is, indeed, in shadow
            continue;
        }

        calcReflection(ray, color, intersectPnt, intersectClr, normal, point_lights[i], incidentRay);
    }

	if(level >= recursionLevel) {
        return tMin;
    }

    return tMin;
}


struct GeneralQuadricSurfaceCoefficient {
    double a, b, c, d, e, f, g, h, i, j;

	GeneralQuadricSurfaceCoefficient(){}

	GeneralQuadricSurfaceCoefficient(double x, double y, double z, double w, double u, double v, double m, double n, double o, double p){
	    a = x;
	    b = y;
	    c = z;
		d = w;
		e = u;
		f = v;
		g = m;
		h = n;
		i = o;
		j = p;
	}

};


class GeneralQuadricSurface: public Object {
    GeneralQuadricSurfaceCoefficient coefficient;
    Vector cubeReferencePoint;
    double length;
    double width;
    double height;

public:
    GeneralQuadricSurface() {
        coefficient.a = 0.0;
		coefficient.b = 0.0;
		coefficient.c = 0.0;
		coefficient.d = 0.0;
		coefficient.e = 0.0;
        coefficient.f = 0.0;
		coefficient.g = 0.0;
		coefficient.h = 0.0;
		coefficient.i = 0.0;
		coefficient.j = 0.0;
        length = 0.0;
		width = 0.0;
	    height = 0.0;
    }

    GeneralQuadricSurface(GeneralQuadricSurfaceCoefficient coeff, Vector cubeRefPoint, double l, double w, double h) {
        coefficient = coeff;
        cubeReferencePoint = cubeRefPoint;
        length = l;
        width = w;
        height = h;
    }

    void draw() {}

    double intersect(Ray, Color&, int);

    ~GeneralQuadricSurface() {}
};

double GeneralQuadricSurface::intersect(Ray ray, Color& color, int level) {
      return -1.0;
}


class Floor: public Object {
    double floorWidth;
    double tileWidth;

    Color foregroundColor;

public:
    Floor() {
        floorWidth = 0.0;
		tileWidth = 0.0;
    }

    Floor(double fWidth, double tWidth, Color fgndColor) {
        floorWidth = fWidth;
        tileWidth = tWidth;
        foregroundColor = fgndColor;
    }

    void draw();
    double intersect(Ray, Color&, int);

    ~Floor() {}
};

void Floor::draw() {

	int row = (int) floorWidth/tileWidth;
	int column = (int) floorWidth/tileWidth;

	for(int i=0; i<row; i++) {
        for(int j=0; j<column; j++) {
            // drawing square on a plane parallel to x-y plane
            if((i+j)%2 == 0)
				glColor3f(getColor().red, getColor().green, getColor().blue);
			else
				glColor3f(foregroundColor.red, foregroundColor.green, foregroundColor.blue);

			double tmpA = -floorWidth/2.0 + tileWidth*j;
			double tmpB = -floorWidth/2.0 + tileWidth*i;
			Vector lbCorner(tmpA, tmpB, 0.0);

            glBegin(GL_QUADS);
            {
                glVertex3f(lbCorner.getX(), lbCorner.getY(), lbCorner.getZ());
                glVertex3f(lbCorner.getX()+tileWidth, lbCorner.getY(), lbCorner.getZ());
                glVertex3f(lbCorner.getX()+tileWidth, lbCorner.getY()+tileWidth, lbCorner.getZ());
                glVertex3f(lbCorner.getX(), lbCorner.getY()+tileWidth, lbCorner.getZ());
            }
            glEnd();
        }
    }
}


double Floor::intersect(Ray ray, Color& color, int level) {

	// determining unit normal vector on appropriate side of floor
    Vector normal(0.0, 0.0, 1.0);

	if(position.dotProduct(normal) > 0.0){
	   normal = normal;
    }
    else{
	   normal = normal.scalarProduct(-1.0);
    }

    // finding intersecting tMin
    double tMin = INF;

    if(normal.dotProduct(ray.getRayDir()) != 0.0) {
        tMin = (-1.0)*(normal.dotProduct(ray.getRayStart()))/(normal.dotProduct(ray.getRayDir()));
    }

	double tmpFW = floorWidth/2.0;

    if(tMin>0.0 && tMin<INF) {

           // make sure the intersection point is on the floor

        Vector intersectPnt = ray.getRayStart() + ray.getRayDir().scalarProduct(tMin);

        if(!((intersectPnt.getX() > -tmpFW && intersectPnt.getX() < tmpFW) && (intersectPnt.getY() > -tmpFW && intersectPnt.getY() < tmpFW))) {
            // intersection point is not on the floor
            tMin = INF;
        }
    }

    if(level == 0) {
        return tMin;
    }

    // illuminating with Phong Lighting Model
    Vector intersectPnt = ray.getRayStart() + ray.getRayDir().scalarProduct(tMin);
    Vector referencePosition = intersectPnt - Vector(-tmpFW, -tmpFW, 0.0);
    Color intersectClr;

	if(((int) (floor(referencePosition.getX()/tileWidth)+floor(referencePosition.getY()/tileWidth)))%2 == 0)
		intersectClr = getColor();
	else
		intersectClr = foregroundColor;

    // calculating ambient light component of reflected ray
    calcAmbient(color, intersectClr);

    // calculating diffuse & specular reflection components of reflected ray
    for(unsigned int i=0; i<point_lights.size(); i++) {
        Ray incidentRay(point_lights[i].getPosition(), intersectPnt-point_lights[i].getPosition());

        // checking if intersection point is in shadow
        double tmp;
		double tMinimum = INF;

        for(unsigned int j=0; j<objects.size(); j++) {
            Color dummyColor;  // black
            tmp = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(tmp>0.0 && tmp<tMinimum) {
                tMinimum = tmp;
            }
        }

        Vector shadowIntersectionPoint = incidentRay.getRayStart() + incidentRay.getRayDir().scalarProduct(tMinimum);
        double epsilon = 0.0000001;  // for tuning light effect

        if(intersectPnt.calcDist(incidentRay.getRayStart())-epsilon > shadowIntersectionPoint.calcDist(incidentRay.getRayStart())) {
            // intersection point is, indeed, in shadow
            continue;
        }

        calcReflection(ray, color, intersectPnt, intersectClr, normal, point_lights[i], incidentRay);
    }

	 if(level >= recursionLevel) {
        return tMin;
    }

    return tMin;
}
