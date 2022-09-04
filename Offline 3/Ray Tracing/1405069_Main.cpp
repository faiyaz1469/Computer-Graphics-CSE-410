#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<limits>

#include<windows.h>
#include<GL/glut.h>

#include "bitmap_image.hpp"
#include "1405069_Header.hpp"

using namespace std;

#define PI 2*acos(0.0)
#define INF numeric_limits<double>::infinity()

double imageWidth;
double imageHeight;

int winWidth = 500;
int winHeight = 500;
double fovY;

bool boolDrawAxes;

int drawaxes;

double rdn_angle = (PI * 3)/180; // 3 degree,in rad format
double _distance = 1;
double axLength = 200;
double viewAngle = (PI * 90)/180;   // 100 degree,in rad format

point camPos;

point u_vec;
point r_vec;
point l_vec;

extern Vector position;
Vector u;
Vector r;
Vector l;

extern int recursionLevel;
int imgPXDim = 0;
int objCnt = 0;
int PointLightsCount = 0;
int SpotLightsCount = 0;

extern vector<Object*> objects;
extern vector<PointLight> point_lights;
extern vector<SpotLight> spot_lights;

int bitmapImgCnt;

/*void drawAxes()
{
	if(drawaxes==1)
	{
	    //glColor3f(1, 1, 1);
		glBegin(GL_LINES);{
		    glColor3f(1, 0, 0);
			glVertex3f( ax_length,0,0);
			glVertex3f(-ax_length,0,0);



			glColor3f(0, 0, 1);
			glVertex3f(0,0, ax_length);
			glVertex3f(0,0,-ax_length);

            glColor3f(0,1,0);
			glVertex3f(0,-ax_length,0);
			glVertex3f(0, ax_length,0);
		}
		glEnd();


	}
}*/

void drawAxes(double axLength) {

	if(!boolDrawAxes) {
        return ;
    }

	glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    {
        glVertex3f(axLength, 0.0, 0.0);
        glVertex3f(-axLength, 0.0, 0.0);
    }
    glEnd();


    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, axLength, 0.0);
        glVertex3f(0.0, -axLength, 0.0);
    }
    glEnd();


    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, 0.0, axLength);
        glVertex3f(0.0, 0.0, -axLength);
    }
    glEnd();
}


void capture() {

    Vector position(camPos.x, camPos.y, camPos.z);
    Vector u(u_vec.x, u_vec.y, u_vec.z);
    Vector r(r_vec.x, r_vec.y, r_vec.z);
    Vector l(l_vec.x, l_vec.y, l_vec.z);

    cout << "Capturing bitmap image at " << "(" << position.getX() << ", " << position.getY() << ", "<< position.getZ() << ")" <<endl;

    bitmap_image bitmapImg(imgPXDim, imgPXDim);

    for(int col=0; col<imgPXDim; col++) {
        for(int row=0; row<imgPXDim; row++) {
            bitmapImg.set_pixel(col, row, 0, 0, 0);  // black
        }
    }

    double planeDistance = winHeight/(2.0*tan(fovY/2.0*PI/180.0));
    Vector topLeft = position + l.scalarProduct(planeDistance) - r.scalarProduct(winWidth/2.0) + u.scalarProduct(winHeight/2.0);

    double du = ((double) winWidth/imgPXDim);
    double dv = ((double) winHeight/imgPXDim);

    topLeft = topLeft + r.scalarProduct(du/2.0) - u.scalarProduct(dv/2.0);

    // capturing the scene from the curr pos of camera

    for(int col=0; col<imgPXDim; col++) {
        for(int row=0; row<imgPXDim; row++) {
            // calculating current pixel and casting ray from camera to (curPixel-camera) direction
            Vector curPixel = topLeft + r.scalarProduct(col*du) - u.scalarProduct(row*dv);
            Ray ray(position, curPixel-position);

            // finding nearest intersecting object (if available) 
            int nearest = INT_MAX;
            double t, tMin=INF;

            for(unsigned int i=0; i<objects.size(); i++) {
                Color color;  // black
                t = objects[i]->intersect(ray, color, 0);

                if(t>0.0 && t<tMin) {
                    tMin = t;
                    nearest = i;
                }
            }

            if(nearest != INT_MAX) {
                Color color;  // black
                tMin = objects[nearest]->intersect(ray, color, 1);
                bitmapImg.set_pixel(col, row, (int) round(color.red*255.0), (int) round(color.green*255.0), (int) round(color.blue*255.0));
            }
        }
    }

    stringstream currBitmapImgCnt;
    currBitmapImgCnt << (++bitmapImgCnt);

    bitmapImg.save_image("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Ray Tracing/"+currBitmapImgCnt.str()+".bmp");
    //cout << "Bitmap image captured at " << "(" <<position.getX() << ", " << position.getY() << ", "<< position.getZ()<<")"<<endl;
	cout << "Bitmap image captured!"<<endl;
}

void RotateOperation(point & fixed, point & rotate1, point & rotate2)
{
    point cross1,cross2;

    cross1 = fixed.cross_product(rotate1);
    cross2 = fixed.cross_product(rotate2);

    rotate1 = rotate1.angle_product(cross1);
    rotate2 = rotate2.angle_product(cross2);

}

void RotateOperation2(point & fixed,point & rotate1,point & rotate2)
{
    point cross1,cross2;

    cross1 = fixed.cross_product(rotate1);
    cross2 = fixed.cross_product(rotate2);

    rotate1 = rotate1.angle_product_2(cross1);
    rotate2 = rotate2.angle_product_2(cross2);

}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

         case '0':
            capture();
            break;

        case '1':
            RotateOperation(u_vec,l_vec,r_vec);
            break;

        case '2':
            RotateOperation2(u_vec,l_vec,r_vec);
            break;

        case '3':
            RotateOperation(r_vec,u_vec,l_vec);
            break;

        case '4':
            RotateOperation2(r_vec,u_vec,l_vec);
            break;

        case '5':
            RotateOperation(l_vec,u_vec,r_vec);
            break;

        case '6':
            RotateOperation2(l_vec,u_vec,r_vec);
            break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			camPos = camPos - l_vec*_distance;
			break;

		case GLUT_KEY_UP:		// up arrow key
            camPos = camPos + l_vec*_distance;
            break;

		case GLUT_KEY_RIGHT:
			camPos = camPos + r_vec*_distance;
            break;

		case GLUT_KEY_LEFT:
			camPos = camPos - r_vec*_distance;
            break;

		case GLUT_KEY_PAGE_UP:
			camPos = camPos + u_vec*_distance;
            break;

		case GLUT_KEY_PAGE_DOWN:
			camPos = camPos - u_vec*_distance;
            break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;

		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y) {
	switch(button) {
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN) {
                boolDrawAxes = !boolDrawAxes;
			}
			break;
        case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN) {

			}
			break;
		default:
			break;
	}
}


void display() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    gluLookAt(camPos.x, camPos.y, camPos.z, camPos.x + l_vec.x, camPos.y + l_vec.y, camPos.z + l_vec.z, u_vec.x, u_vec.y, u_vec.z);
	glMatrixMode(GL_MODELVIEW);

	drawAxes(300.0);

	for(unsigned int i=0; i<objects.size(); i++) {
        objects[i]->draw();
	}

	for(unsigned int i=0; i<point_lights.size(); i++) {
        point_lights[i].draw();
	}

	for(unsigned int i=0; i<spot_lights.size(); i++) {
        spot_lights[i].draw();
	}

	glutSwapBuffers();
}

void animate() {
	glutPostRedisplay();
}

void init() {

	boolDrawAxes = true;
	drawaxes=1;

	fovY = 80.0;

    camPos.x = 180;
	camPos.y = 140;
	camPos.z = 50;

	u_vec.x = 0;
	u_vec.y = 0;
	u_vec.z = 1;

	r_vec.x = -1/sqrt(2);
	r_vec.y = 1/sqrt(2);
	r_vec.z = 0;

	l_vec.x = -1/sqrt(2);
	l_vec.y = -1/sqrt(2);
	l_vec.z = 0;

	bitmapImgCnt = 0;

	glClearColor(0, 0, 0, 0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovY, 1.0, 1.0, 1000.0);
}


void freeObjects() {
    for(unsigned int i=0; i<objects.size(); i++) {
        delete objects[i];
    }

    objects.clear();
    point_lights.clear();
	spot_lights.clear();
	cout<<"Memory freed!"<<endl;
}

void loadData() {

    ifstream input("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Ray Tracing/scene.txt",ios::in);

    input >> recursionLevel >> imgPXDim;
    input >> objCnt;

    string objectShape;
    Object* object = NULL;

    for(int i=0; i<objCnt; i++) {
        input >> objectShape;

        if(objectShape == "sphere"){
            double radius;
            double centerX,centerY,centerZ;

            input >> centerX >> centerY >> centerZ;
            input >> radius;

            Vector center(centerX,centerY,centerZ);

            object = new Sphere(center, radius, 72, 24);
        }
        else if(objectShape == "triangle") {

           double aX,aY,aZ,bX,bY,bZ,cX,cY,cZ;

            input >> aX >> aY >> aZ;
            input >> bX >> bY >> bZ;
            input >> cX >> cY >> cZ;

            Vector a(aX,aY,aZ);
            Vector b(bX,bY,bZ);
            Vector c(cX,cY,cZ);

            object = new Triangle(a, b, c);
        }
        else if(objectShape == "general") {

            double coA,coB,coC,coD,coE,coF,coG,coH,coI,coJ;
            double refX,refY,refZ;
            double length, width, height;

            input >> coA >> coB >> coC >> coD >> coE >> coF >> coG >> coH >> coI >> coJ;
            input >> refX >> refY >> refZ;
            input >> length >> width >> height;

            GeneralQuadricSurfaceCoefficient coefficient(coA,coB,coC,coD,coE,coF,coG,coH,coI,coJ);
            Vector cubeReferencePoint(refX,refY,refZ);

            object = new GeneralQuadricSurface(coefficient, cubeReferencePoint, length, width, height);
        }

        double r,g,b;
        double amb,diff,spec,recur;
        int shininess;

        input >> r >> g >> b;
        input >> amb >> diff >> spec >> recur;
        input >> shininess;

        Color color(r,g,b);
        ReflectionCoefficient reflectionCoefficient(amb,diff,spec,recur);

        object->setColor(color);
        object->setReflectionCoefficient(reflectionCoefficient);
        object->setShininess(shininess);

        objects.push_back(object);
    }
    object = NULL;

    input >> PointLightsCount;

    for(int i=0; i<PointLightsCount; i++) {

        double posX,posY,posZ;
        double r,g,b;

        input >> posX >> posY >> posZ;
		//cout<<posX <<" " << posY<< " " <<posZ<<endl;
        input >> r >> g >> b;
		//cout<<r <<" " << g<< " " <<b<<endl;

         Vector position(posX,posY,posZ);
         Color color(r,g,b);

        point_lights.push_back(PointLight(position, color, 1.0, 12, 4));
    }

	input >> SpotLightsCount;

    for(int i=0; i<SpotLightsCount; i++) {

        double posX,posY,posZ;
        double r,g,b;
		double dirX,dirY,dirZ;
		double coAngle;

        input >> posX >> posY >> posZ;
		//cout<<posX <<" " << posY<< " " <<posZ<<endl;
        input >> r >> g >> b;
		//cout<<r <<" " << g<< " " <<b<<endl;
		input >> dirX >> dirY >> dirZ;
		//cout<<dirX <<" " << dirY<< " " <<dirZ<<endl;
		input >> coAngle;
		//cout<< coAngle <<endl;

         Vector position(posX,posY,posZ);
         Color color(r,g,b);
		 PointLight point_light(position,color);
		 Vector direction(dirX, dirY, dirZ);

		spot_lights.push_back(SpotLight(point_light, direction, coAngle, 1.0, 12, 4));
    }
    input.close();

	cout<<"File Reading completed!"<<endl;

    object = new Floor(1000.0, 20.0, Color());  // black

    object->setColor(Color(1.0, 1.0, 1.0));  // white
    object->setReflectionCoefficient(ReflectionCoefficient(0.25, 0.25, 0.25, 0.25));
    object->setShininess(15);

    objects.push_back(object);
    object = NULL;
}


int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitWindowSize(winWidth, winHeight);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

	glutCreateWindow("Ray Casting & Illumination");

	init();

	glEnable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
	glutIdleFunc(animate);

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);


	loadData();
    glutMainLoop();
    freeObjects();

	return 0;
}
