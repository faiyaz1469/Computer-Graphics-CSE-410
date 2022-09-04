#include <bits/stdc++.h>
#include<windows.h>
#include<GL/glut.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

using namespace std;

struct point
{
	double x,y,z;
    point(){
	   x=0;
	   y=0;
	   z=0;
	}

	point(double x1, double y1, double z1){
        x = x1;
        y = y1;
        z = z1;
    }
};

point camPos;
point u_vec;
point r_vec;
point l_vec;

double camHeight;
double camAngle;
int drawgrid ;
int drawaxes;
double angle;

double _z_angle = 0;
double _x_angle = 0;
double _x_cdr_angle = 0;
double _own_angle = 0;
double rdn_angle = (pi * 3)/180;

double _distance = 3;
double width = 35;
double rad = 15;

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);

		}glEnd();
	}
}

void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f(-a, a,0);
		glVertex3f(-a,-a,0);
		glVertex3f( a,-a,0);
    }glEnd();
}

void drawSphere(double rad,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=rad*sin(((double)i/(double)stacks)*(pi/2));
		r=rad*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*pi/2);
			points[i][j].y=r*sin(((double)j/(double)slices)*pi/2);
			points[i][j].z=h;
		}
	}
	glColor3f(1.0,0,0);
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawCylinder(double rad,double height,int segments)
{
    int i;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=rad*cos(((double)i/(double)segments)*pi/2);
        points[i].y=rad*sin(((double)i/(double)segments)*pi/2);
    }
    glColor3f(0,1.0,0);
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {

        glBegin(GL_QUADS);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,height);
			glVertex3f(points[i].x,points[i].y,height);
        }
        glEnd();
    }
}

void drawUpperCorner(double w, double r){
    glPushMatrix();
        glTranslatef(w-r,w-r,w-r);
        drawSphere(r,20,20);
    glPopMatrix();
}

void drawLowerCorner(double w, double r){
    glPushMatrix();
        glTranslatef(w-r,w-r,-(w-r));
        glRotatef(90, 0,1,0);
        drawSphere(r,20,20);
    glPopMatrix();
}

void drawHorizontalEdge(double w,double r){
    glPushMatrix();
        glRotatef(-90, 0,1,0);
        glTranslatef(w-r,w-r,-(w-r));
        drawCylinder(r,2*(w-r),20);
    glPopMatrix();
}

void drawVerticalEdge(double w,double r){
    glPushMatrix();
        glTranslatef(w-r,w-r,-(w-r));
        drawCylinder(r,2*(w-r),20);
    glPopMatrix();
}

void drawSurface(double a){
    a = abs(a);
    glColor3f(1.0,1.0,1.0);
    drawSquare(a);
}

void drawSS()
{
    //frontCorner
    drawVerticalEdge(width,rad);

    //backCorner
    glPushMatrix();
        glRotatef(180, 0,0,1);
        drawVerticalEdge(width,rad);
    glPopMatrix();

    //rightCorner
    glPushMatrix();
        glRotatef(90, 0,0,1);
        drawVerticalEdge(width,rad);
    glPopMatrix();

    //leftCorner
    glPushMatrix();
        glRotatef(-90, 0,0,1);
        drawVerticalEdge(width,rad);
    glPopMatrix();

    drawHorizontalEdge(width,rad);     //upperFrontRight

    glPushMatrix();
        glRotatef(90, 1,0,0);
        drawHorizontalEdge(width,rad);  //upperBackRight
    glPopMatrix();

    glPushMatrix();
        glRotatef(180, 1,0,0);
        drawHorizontalEdge(width,rad);    //lowerBackRight
    glPopMatrix();

    glPushMatrix();
        glRotatef(-90, 1,0,0);
        drawHorizontalEdge(width,rad);   //lowerFrontRight
    glPopMatrix();

    glPushMatrix();
        glRotatef(90, 0,0,1);
        drawHorizontalEdge(width,rad);    //upperBackLeft
    glPopMatrix();

    glPushMatrix();
        glRotatef(90, 0,1,0);
        glRotatef(90, 0,0,1);
        drawHorizontalEdge(width,rad);   //upperFrontLeft
    glPopMatrix();

    glPushMatrix();
        glRotatef(180, 0,1,0);
        glRotatef(90, 0,0,1);
        drawHorizontalEdge(width,rad);    //lowerFrontLeft
    glPopMatrix();

    glPushMatrix();
        glRotatef(-90, 0,1,0);
        glRotatef(90, 0,0,1);
        drawHorizontalEdge(width,rad);     //lowerBackLeft
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0,0,width);
        drawSurface(width-rad);     //Top
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0,0,-width);
        drawSurface(width-rad);     //Bottom
    glPopMatrix();


    glColor3f(1,0,0);
    glPushMatrix();
        glTranslatef(width,0,0);
        glRotatef(90, 0,1,0);
        drawSurface(width-rad);     //frontLeft
    glPopMatrix();

    glPushMatrix();
        glTranslatef(-width,0,0);
        glRotatef(90, 0,1,0);
        drawSurface(width-rad);    //backLeft
    glPopMatrix();

    glColor3f(0,1,0);
    glPushMatrix();
        glTranslatef(0,width,0);
        glRotatef(90, 1,0,0);
        drawSurface(width-rad);   //frontRight
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0,-width,0);
        glRotatef(90, 1,0,0);
        drawSurface(width-rad);    //backRight
    glPopMatrix();

    //frontCorner
    drawUpperCorner(width,rad);
    drawLowerCorner(width,rad);

    //backCorner
    glPushMatrix();
        glRotatef(180, 0,0,1);
        drawUpperCorner(width,rad);
        drawLowerCorner(width,rad);
    glPopMatrix();

    //rightCorner
    glPushMatrix();
        glRotatef(90, 0,0,1);
        drawUpperCorner(width,rad);
        drawLowerCorner(width,rad);
    glPopMatrix();

    //leftCorner
    glPushMatrix();
        glRotatef(-90, 0,0,1);
        drawUpperCorner(width,rad);
        drawLowerCorner(width,rad);
    glPopMatrix();
}


void op_rotate(point & fx, point & rt1, point & rt2, double _rdn)
{

    point cross1,cross2;

    cross1.x = fx.y * rt1.z - fx.z * rt1.y;
    cross1.y = fx.z * rt1.x - fx.x * rt1.z;
    cross1.z = fx.x * rt1.y - fx.y * rt1.x;

    cross2.x = fx.y * rt2.z - fx.z * rt2.y;
    cross2.y = fx.z * rt2.x - fx.x * rt2.z;
    cross2.z = fx.x * rt2.y - fx.y * rt2.x;

    rt1.x = rt1.x * cos(_rdn) + cross1.x * sin(_rdn);
    rt1.y = rt1.y * cos(_rdn) + cross1.y * sin(_rdn);
    rt1.z = rt1.z * cos(_rdn) + cross1.z * sin(_rdn);

    rt2.x = rt2.x * cos(_rdn) + cross2.x * sin(_rdn);
    rt2.y = rt2.y * cos(_rdn) + cross2.y * sin(_rdn);
    rt2.z = rt2.z * cos(_rdn) + cross2.z * sin(_rdn);

}

void keyboardListener(unsigned char key, int x,int y){

	switch(key){
         case '1':
            op_rotate(u_vec,l_vec,r_vec,rdn_angle);
            break;

        case '2':
            op_rotate(u_vec,l_vec,r_vec,-rdn_angle);
            break;

        case '3':
            op_rotate(r_vec,u_vec,l_vec,rdn_angle);
            break;
        case '4':
            op_rotate(r_vec,u_vec,l_vec,-rdn_angle);
            break;
        case '5':
            op_rotate(l_vec,u_vec,r_vec,rdn_angle);
            break;
        case '6':
            op_rotate(l_vec,u_vec,r_vec,-rdn_angle);
            break;

		default:
			break;
	}
}

void specialKeyListener(int key, int x,int y){

	switch(key){

		case GLUT_KEY_DOWN:
			camPos.x -= _distance*l_vec.x;
			camPos.y -= _distance*l_vec.y;
			camPos.z -= _distance*l_vec.z;

            break;

		case GLUT_KEY_UP:
			camPos.x += _distance*l_vec.x;
			camPos.y += _distance*l_vec.y;
			camPos.z += _distance*l_vec.z;

			break;

		case GLUT_KEY_RIGHT:
			camPos.x += _distance*r_vec.x;
			camPos.y += _distance*r_vec.y;
			camPos.z += _distance*r_vec.z;

			break;

		case GLUT_KEY_LEFT:
			camPos.x -= _distance*r_vec.x;
			camPos.y -= _distance*r_vec.y;
			camPos.z -= _distance*r_vec.z;

			break;

		case GLUT_KEY_PAGE_UP:
			camPos.x += _distance*u_vec.x;
			camPos.y += _distance*u_vec.y;
			camPos.z += _distance*u_vec.z;

			break;

		case GLUT_KEY_PAGE_DOWN:
			camPos.x -= _distance*u_vec.x;
			camPos.y -= _distance*u_vec.y;
			camPos.z -= _distance*u_vec.z;

			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
            rad += 0.8;
            if(rad >= width)
                rad = width;
			break;

		case GLUT_KEY_END:
		    rad -= 0.8;
            if(rad <= 0.8)
                rad = 0.0;
			break;

		default:
			break;
	}
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)

	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes = 1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
    gluLookAt(camPos.x, camPos.y, camPos.z, camPos.x + l_vec.x, camPos.y + l_vec.y, camPos.z + l_vec.z, u_vec.x, u_vec.y, u_vec.z);

	glMatrixMode(GL_MODELVIEW);
	drawAxes();
	drawGrid();
    drawSS();
	glutSwapBuffers();
}

void animate(){
	angle+=0.05;
	glutPostRedisplay();
}

void init(){

	drawgrid=0;
	drawaxes=1;
	camHeight=150.0;
	camAngle=1.0;
	angle=0;

	camPos = point(100, 100, 0);
	u_vec = point(0, 0, 1);
	r_vec = point(-1/sqrt(2), 1/sqrt(2), 0);
	l_vec = point(-1/sqrt(2), -1/sqrt(2), 0);

	glClearColor(0, 0, 0, 0);
    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(80, 1, 1, 1000.0);
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My Offline Sphere To/From Cube with Fully Controllable Camera");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
