#include <bits/stdc++.h>
#include<windows.h>
#include<GL/glut.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<iostream>

#define pi (2*acos(0.0))

double camHeight;
double camAngle;
int drawgrid;
int drawaxes;
double angle;
double dist = 3.0;
double shiftAngle = 0.03;

double moveDist = 1;
double rdn_angle = (pi * 2)/180;

double wheelX,wheelA,wheelRad,wheelDir;

struct point{
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

point pos;
point dir;
point fixed;

void move1(double);
void rotate1(double);

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
		glColor3f(0.4, 0.4, 0.4);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				//if(i==0)
				//	continue;	//SKIP the MAIN axes

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

void drawWheel(double radius,double width){

    int segments = 32;
    int i;
    double shade;
    struct point points[100];

    //generate points
    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z + wheelRad);
    glRotatef(atan2(-dir.y,-dir.x)/pi*180, 0, 0, 1);
    glRotatef(90, 1,0,0);
    glRotatef(wheelA, 0,0,1);

    for(i=0;i<=segments;i++)
    {
        points[i].x = radius*cos(((double)i/(double)segments)*pi*2);
        points[i].y = radius*sin(((double)i/(double)segments)*pi*2);
    }

    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2){
            shade=2*(double)i/(double)segments;
        }
        else{
            shade=2*(1.0-(double)i/(double)segments);
        }

        glColor3f(shade,shade,shade);

        glBegin(GL_QUADS);
        {
			glVertex3f(points[i].x, points[i].y, -width/1.5);
			glVertex3f(points[i+1].x, points[i+1].y, -width/1.5);
			glVertex3f(points[i+1].x, points[i+1].y, width/1.5);
			glVertex3f(points[i].x, points[i].y, width/1.5);
        }
        glEnd();
    }

    glBegin(GL_QUADS);
    {
        glColor3f(.6,.6,.6);
        glVertex3f(radius,0,-width/4);
        glVertex3f(-radius,0,-width/4);
        glVertex3f(-radius,0, width/4);
        glVertex3f(radius,0, width/4);

        glVertex3f(0,radius,-width/4);
        glVertex3f(0,-radius,-width/4);
        glVertex3f(0,-radius, width/4);
        glVertex3f(0, radius, width/4);
    }
    glEnd();
    glPopMatrix();
}


point op_rotate(point & l, point & r, double _rdn)
{
    point u,ans ;

    u.x = r.y*l.z - l.y*r.z;
    u.y = r.z*l.x - l.z*r.x;
    u.z = r.x*l.y - l.x*r.y;

    ans.x = l.x * cos(_rdn) + u.x * sin(_rdn);
    ans.y = l.y * cos(_rdn) + u.y * sin(_rdn);
    ans.z = l.z * cos(_rdn) + u.z * sin(_rdn);

    return ans;

}

void move1(double slide){

    pos.x = pos.x + dir.x*slide;
    pos.y = pos.y + dir.y*slide;
    pos.z = pos.z + dir.z*slide;

    double dA = 360.0/(2*pi*wheelRad)*slide;
    wheelA += dA;

    if(wheelA > 360)
        wheelA -= 360;
}

void rotate1(double angle){
    dir = op_rotate(dir, fixed, angle );
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

        case 'w':
			move1(moveDist);
			break;

        case 's':
			move1(-moveDist);
			break;

        case 'a':
			rotate1(rdn_angle);
			break;

        case 'd':
			rotate1(-rdn_angle);
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){

		case GLUT_KEY_DOWN:
			camHeight -= dist;
			break;

		case GLUT_KEY_UP:
			camHeight += dist;
			break;

		case GLUT_KEY_RIGHT:
			camAngle += shiftAngle;
			break;

		case GLUT_KEY_LEFT:
			camAngle -= shiftAngle;
			break;

		case GLUT_KEY_PAGE_UP:
			break;

		case GLUT_KEY_PAGE_DOWN:
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


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			//........
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

	gluLookAt(80*cos(camAngle), 80*sin(camAngle), camHeight,		0,0,0,		0,0,1);

	glMatrixMode(GL_MODELVIEW);

    drawAxes();
	drawGrid();

	drawWheel(wheelRad,wheelRad/3);
	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){

	drawgrid=1;
	drawaxes=0;
	camHeight=60.0;
	camAngle=1.0;
	angle=0;

	glClearColor(0,0,0,0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,	1,	1,	1000.0);

    wheelRad = 10.5;
	wheelX = wheelDir = wheelA = 0;

	pos = point( 0, 0, 0);
    dir = point(-1, 0, 0);
    fixed = point(0, 0, 1);
}

int main(int argc, char **argv){

	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My Offline Wheel");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
