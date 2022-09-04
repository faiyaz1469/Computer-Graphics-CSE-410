#include <bits/stdc++.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<string>
#include<cstdlib>
#include<stack>
#include<ctime>
#include<limits>
#include "bitmap_image.hpp"

using namespace std;

#define PI 2.0*acos(0.0)
#define INF numeric_limits<double>::infinity()

double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, near, far;
int triagleCount=0, pushCounter=0;
int screenWidth, screenHeight;
double leftLimitX, rightLimitX, bottomLimitY, topLimitY, frontLimitZ, rearLimitZ;
double dx, dy, topY, bottomY, leftX, rightX;


class Point {

    double x;
    double y;
    double z;
    double w;

public:

    Point() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 1.0;
    }

    Point(double x1, double y1, double z1) {
        x = x1;
        y = y1;
        z = z1;
        w = 1;
    }

    Point(double x1, double y1, double z1, double w1) {
        x = x1;
        y = y1;
        z = z1;
        w = w1;
    }

    double getX() const {
        return x;
    }

    void setX(double x1) {
        x = x1;
    }

    double getY() const {
        return y;
    }

    void setY(double y1) {
        y = y1;
    }

    double getZ() const {
        return z;
    }

    void setZ(double z1) {
        z = z1;
    }

    double getW() const {
        return w;
    }

    void normalize()
    {
        double r = sqrt((x * x) + (y * y) + (z * z));
        x = x / r;
        y = y / r;
        z = z / r;
    }

    void scale()
    {
        x  = x / w;
        y  = y / w;
        z  = z / w;
        w  = w / w;
    }

    friend ifstream& operator>>(ifstream&, Point&);
    friend ofstream& operator<<(ofstream&, Point&);

    ~Point() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 1.0;
    }
};

ifstream& operator>>(ifstream &input, Point &point) {
    input >> point.x >> point.y >> point.z;
    return input;
}

ofstream& operator<<(ofstream &output, Point &point) {
    output << fixed << setprecision(7) << point.x << ' ' << point.y << ' '  << point.z;
    return output;
}


struct Color {
    int redVal;
    int greenVal;
    int blueVal;
};

struct Triangle {
    Point corners[3];
    Color rgb;
};


class Transformation {

    double mat[4][4];

    void createIdentityMatrix(){

        int matrixrow = 4;
        int matrixcol = 4;

        for(int i=0; i<matrixrow; i++) {
            for(int j=0; j<matrixcol; j++) {
                if(i==j){
                    mat[i][j] = 1.0;
                }
                else{
                    mat[i][j] = 0.0;
                }
            }
        }
    }

    Point applyRodriguesFormula (Point x, Point a, double radian) {

        double dot = x.getX() * a.getX() + x.getY() * a.getY() + x.getZ() * a.getZ();

        double crossX = a.getY() * x.getZ() - a.getZ() * x.getY();
        double crossY = a.getZ() * x.getX() - a.getX() * x.getZ();
        double crossZ = a.getX() * x.getY() - a.getY() * x.getX();

        Point cross(crossX,crossY,crossZ);

        double productX = x.getX()*cos(radian) + a.getX()*(dot)*(1-cos(radian)) + cross.getX()*sin(radian);
        double productY = x.getY()*cos(radian) + a.getY()*(dot)*(1-cos(radian)) + cross.getY()*sin(radian);
        double productZ = x.getZ()*cos(radian) + a.getZ()*(dot)*(1-cos(radian)) + cross.getZ()*sin(radian);

        Point product(productX,productY,productZ);

        return product;
    }


public:

    Transformation() {
        createIdentityMatrix();
    }

    void multMatrix (Transformation transformation){

        double temp[4][4];

        for(int i=0; i<4; i++) {
            for(int j=0; j<4; j++) {
                temp[i][j] = 0.0;

                for(int k=0; k<4; k++) {
                    temp[i][j] += mat[i][k]*transformation.mat[k][j];
                }
            }
        }
        for(int i=0; i<4; i++) {
            for(int j=0; j<4; j++) {
                mat[i][j] = temp[i][j];
            }
        }
    }

    void createTranslationMatrix (double tx, double ty, double tz){

        mat[0][3] = tx;
        mat[1][3] = ty;
        mat[2][3] = tz;
    }

    void createScalingMatrix (double sx, double sy, double sz){

        mat[0][0] = sx;
        mat[1][1] = sy;
        mat[2][2] = sz;
    }

    void createRotationMatrix(double, double, double, double);
    void createViewMatrix(Point, Point, Point);
    void createProjectionMatrix(double, double, double, double);

    Point operator*(const Point);

    ~Transformation() {
        createIdentityMatrix();
    }
};


void Transformation::createRotationMatrix(double angle, double ax, double ay, double az) {

    Point a(ax, ay, az);
    Point i(1.0, 0.0, 0.0);
    Point j(0.0, 1.0, 0.0);
    Point k(0.0, 0.0, 1.0);

    a.normalize();

    double rad = (PI * angle)/180;
    Point c1 = applyRodriguesFormula(i, a, rad);
    Point c2 = applyRodriguesFormula(j, a, rad);
    Point c3 = applyRodriguesFormula(k, a, rad);

    mat[0][0] = c1.getX();
    mat[1][0] = c1.getY();
    mat[2][0] = c1.getZ();

    mat[0][1] = c2.getX();
    mat[1][1] = c2.getY();
    mat[2][1] = c2.getZ();

    mat[0][2] = c3.getX();
    mat[1][2] = c3.getY();
    mat[2][2] = c3.getZ();
}


void Transformation::createViewMatrix(Point eye, Point look, Point up) {

    double lX = look.getX()-eye.getX();
    double lY = look.getY()-eye.getY();
    double lZ = look.getZ()-eye.getZ();

    Point l(lX,lY,lZ);
    l.normalize();

    double rX = l.getY() * up.getZ() - l.getZ() * up.getY();
    double rY = l.getZ() * up.getX() - l.getX() * up.getZ();
    double rZ = l.getX() * up.getY() - l.getY() * up.getX();

    Point r(rX,rY,rZ);
    r.normalize();

    double uX = r.getY() * l.getZ() - r.getZ() * l.getY();
    double uY = r.getZ() * l.getX() - r.getX() * l.getZ();
    double uZ = r.getX() * l.getY() - r.getY() * l.getX();

    Point u(uX,uY,uZ);

    Transformation trans_Transform;
    trans_Transform.createTranslationMatrix(-eye.getX(), -eye.getY(), -eye.getZ());

    mat[0][0] = r.getX();
    mat[0][1] = r.getY();
    mat[0][2] = r.getZ();

    mat[1][0] = u.getX();
    mat[1][1] = u.getY();
    mat[1][2] = u.getZ();

    mat[2][0] = -l.getX();
    mat[2][1] = -l.getY();
    mat[2][2] = -l.getZ();

    multMatrix(trans_Transform);
}


void Transformation::createProjectionMatrix(double fovY, double aspectRatio, double near, double far) {

    double fovX = fovY*aspectRatio;
    double rad1 = ((fovX/2)*PI)/180;
    double rad2 = ((fovY/2)*PI)/180;

    double r = near*tan(rad1);
    double t = near*tan(rad2);

    mat[0][0] = near/r;
    mat[1][1] = near/t;
    mat[2][2] = -(far+near)/(far-near);
    mat[2][3] = -(2.0*far*near)/(far-near);
    mat[3][2] = -1.0;
    mat[3][3] = 0.0;
}


Point Transformation::operator*(const Point point) {

            double temp[4];

            for(int i=0; i<4; i++) {
                temp[i] = 0.0;

                for(int j=0; j<4; j++) {

                        if(j==0)
                           temp[i] += mat[i][j]*point.getX();
                        else if(j==1)
                           temp[i] += mat[i][j]*point.getY();
                        else if(j==2)
                           temp[i] += mat[i][j]*point.getZ();
                        else
                           temp[i] += mat[i][j]*point.getW();
                }
            }

            double x = temp[0];
            double y = temp[1];
            double z = temp[2];
            double w = temp[3];

            Point result(x, y, z, w);

            return result;
}


void readScene(ifstream & input)
{
    input >> eyeX >> eyeY >> eyeZ;
    input >> lookX >> lookY >> lookZ;
    input >> upX >> upY >> upZ;
    input >> fovY >> aspectRatio >> near >> far;
}


void modelTransForm(ifstream & input, ofstream & output)
{
    string cmd;

    stack<Transformation> transform_MatStack;
    transform_MatStack.push(Transformation());

    while(true) {

        input >> cmd;

        if(cmd.compare("triangle") == 0) {

            Point p1, p2, p3;

            input >> p1;
            input >> p2;
            input >> p3;

            p1 = transform_MatStack.top()*p1;
            p1.scale();

            p2 = transform_MatStack.top()*p2;
            p2.scale();

            p3 = transform_MatStack.top()*p3;
            p3.scale();

            output << p1 << endl;
            output << p2 << endl;
            output << p3 << endl;
            output << endl;

            triagleCount++;

        }
        else if(cmd.compare("translate") == 0) {

            double tx, ty, tz;
            input >> tx >> ty >> tz;

            Transformation translationTransformation;
            translationTransformation.createTranslationMatrix(tx, ty, tz);

            Transformation temp;
            Transformation top_Transform = transform_MatStack.top();

            temp.multMatrix(top_Transform);
            temp.multMatrix(translationTransformation);

            transform_MatStack.pop();
            transform_MatStack.push(temp);

        }
        else if(cmd.compare("scale") == 0) {

            double sx, sy, sz;
            input >> sx >> sy >> sz;

            Transformation scalingTransformation;
            scalingTransformation.createScalingMatrix(sx, sy, sz);

            Transformation temp;
            Transformation top_Transform = transform_MatStack.top();

            temp.multMatrix(top_Transform);
            temp.multMatrix(scalingTransformation);

            transform_MatStack.pop();
            transform_MatStack.push(temp);

        }
        else if(cmd.compare("rotate") == 0) {

            double angle, ax, ay, az;
            input >> angle >> ax >> ay >> az;

            Transformation rotationTransformation;
            rotationTransformation.createRotationMatrix(angle, ax, ay, az);

            Transformation temp;
            Transformation topTransformation = transform_MatStack.top();

            temp.multMatrix(topTransformation);
            temp.multMatrix(rotationTransformation);

            transform_MatStack.pop();
            transform_MatStack.push(temp);

        }
        else if(cmd.compare("push") == 0) {

            Transformation temp;
            Transformation topTransformation = transform_MatStack.top();

            temp.multMatrix(temp);
            temp.multMatrix(topTransformation);

            transform_MatStack.push(temp);
            pushCounter++;

        }
        else if(cmd.compare("pop") == 0) {

            if(pushCounter == 0) {
                cout << cmd << ": pop on empty stack" << endl;
                break;
            }

            transform_MatStack.pop();
            pushCounter--;

        }
        else if(cmd.compare("end") == 0) {

            break;

        }
    }

    cout <<"Stage1 completed!"<<endl;

    input.close();
    output.close();
}


void viewTransForm(ifstream & input, ofstream & output)
{
    Point e(eyeX, eyeY, eyeZ);
    Point l(lookX, lookY, lookZ);
    Point u(upX, upY, upZ);

    Transformation view_Transform;
    view_Transform.createViewMatrix(e,l,u);

    for(int i=0; i<triagleCount; i++) {
        Point p1, p2, p3;

        input >> p1;
        input >> p2;
        input >> p3;

        p1 = view_Transform*p1;
        p1.scale();

        p2 = view_Transform*p2;
        p2.scale();

        p3 = view_Transform*p3;
        p3.scale();

        output << p1 << endl;
        output << p2 << endl;
        output << p3 << endl;
        output << endl;
    }

    cout <<"Stage2 completed!"<<endl;

    input.close();
    output.close();
}


void projectionTransForm(ifstream & input, ofstream & output)
{

    Transformation project_Transform;
    project_Transform.createProjectionMatrix(fovY, aspectRatio, near, far);

    for(int i=0; i<triagleCount; i++) {
        Point p1, p2, p3;

        input >> p1;
        input >> p2;
        input >> p3;

        p1 = project_Transform*p1;
        p1.scale();

        p2 = project_Transform*p2;
        p2.scale();

        p3 = project_Transform*p3;
        p3.scale();

        output << p1 << endl;
        output << p2 << endl;
        output << p3 << endl;
        output << endl;
    }

    cout <<"Stage3 completed!"<<endl;

    input.close();
    output.close();
}


vector <Triangle> tri_angle;
double** z_Buff;
Color** frame_Buff;

void readConfig(ifstream & input)
{

    input >> screenWidth >> screenHeight;
    input >> leftLimitX;
    input >> bottomLimitY;
    input >> frontLimitZ >> rearLimitZ;

    input.close();

    rightLimitX = -leftLimitX;
    topLimitY = -bottomLimitY;
}

void readStage3(ifstream & input)
{
    srand(time(0));

    for(int i = 0; i < triagleCount; i++) {
        input >> tri_angle[i].corners[0];
        input >> tri_angle[i].corners[1];
        input >> tri_angle[i].corners[2];

        tri_angle[i].rgb.redVal = rand()%256;
        tri_angle[i].rgb.greenVal = rand()%256;
        tri_angle[i].rgb.blueVal = rand()%256;
    }
    input.close();
}


void initZ_buffer(){

    dx = (rightLimitX - leftLimitX)/screenWidth;
    dy = (topLimitY - bottomLimitY)/screenHeight;

    topY = topLimitY - (dy/2.0);
    leftX = leftLimitX + (dx/2.0);

    bottomY = bottomLimitY + (dy/2.0);
    rightX = rightLimitX - (dx/2.0);

    z_Buff = new double*[screenHeight];
    for(int i=0; i<screenHeight; i++) {
        z_Buff[i] = new double[screenWidth];
    }

    for(int row=0; row<screenHeight; row++) {
        for(int col=0; col<screenWidth; col++) {
            z_Buff[row][col] = rearLimitZ;
        }
    }

    frame_Buff = new Color*[screenHeight];
    for(int i=0; i<screenHeight; i++) {
        frame_Buff[i] = new Color[screenWidth];
    }

    for(int row=0; row<screenHeight; row++) {
        for(int col=0; col<screenWidth; col++) {
            frame_Buff[row][col].redVal = 0;
            frame_Buff[row][col].greenVal = 0;
            frame_Buff[row][col].blueVal = 0;
        }
    }
}


void applyProcedure(){

    for(int i=0, topScan_Line, bottomScan_Line; i<triagleCount; i++) {
        // finding topScan_Line & bottomScan_Line after clipping
        double maxY = max(tri_angle[i].corners[0].getY(), max(tri_angle[i].corners[1].getY(), tri_angle[i].corners[2].getY()));
        double minY = min(tri_angle[i].corners[0].getY(), min(tri_angle[i].corners[1].getY(), tri_angle[i].corners[2].getY()));

        if(maxY >= topY) {
            topScan_Line = 0;
        }
        else {
            topScan_Line = (int) round((topY - maxY)/dy);
        }

        if(minY <= bottomY) {
            bottomScan_Line = screenHeight - 1;
        }
        else {
            bottomScan_Line = screenHeight - (1 + ((int) round((minY - bottomY)/dy)));
        }

        // scanning from topScan_Line to bottomScan_Line (incl)
        for(int row = topScan_Line, leftIntersectCol, rightIntersectCol; row <= bottomScan_Line; row++) {

            double ys = topY - row*dy;

            Point intersectPts[3];
            intersectPts[0] = Point(INF, ys, 0, 1);
            intersectPts[1] = Point(INF, ys, 1, 2);
            intersectPts[2] = Point(INF, ys, 2, 0);

           for(int j=0; j<3; j++) {

                Point p1 = tri_angle[i].corners[(int) intersectPts[j].getZ()];
                Point p2 = tri_angle[i].corners[(int) intersectPts[j].getW()];

                if(p1.getY() != p2.getY()) {
                    intersectPts[j].setX(p1.getX() + (ys - p1.getY())*(p1.getX() - p2.getX())/(p1.getY() - p2.getY()));
                }
            }

            // filtering out all invalid points
            for(int j=0; j<3; j++) {

                Point p1 = tri_angle[i].corners[(int) intersectPts[j].getZ()];
                Point p2 = tri_angle[i].corners[(int) intersectPts[j].getW()];

                if(intersectPts[j].getX() != INF) {
                    if(intersectPts[j].getX() > max(p1.getX(), p2.getX()) || intersectPts[j].getX() < min(p1.getX(), p2.getX()) || intersectPts[j].getY() > max(p1.getY(), p2.getY()) || intersectPts[j].getY() < min(p1.getY(), p2.getY())) {
                        intersectPts[j].setX(INF);
                    }
                }
            }

            // finding out leftIntersecting & rightIntersecting points
            int maxIndex = -1, minIndex = -1;
            double maxX, minX;

            for(int j=0; j<3; j++) {

                if(maxIndex == -1 && minIndex == -1) {
                    if(intersectPts[j].getX() != INF) {
                        maxIndex = minIndex = j;
                        maxX = minX = intersectPts[j].getX();
                    }
                }
                else {
                    if(intersectPts[j].getX() != INF) {
                        if(intersectPts[j].getX() < minX) {
                            minIndex = j;
                            minX = intersectPts[j].getX();
                        }
                        if(intersectPts[j].getX() > maxX) {
                            maxIndex = j;
                            maxX = intersectPts[j].getX();
                        }
                    }
                }
            }

            // finding leftIntersectCol & rightIntersectCol after clipping
            if(intersectPts[minIndex].getX() <= leftX) {
                leftIntersectCol = 0;
            }
            else {
                leftIntersectCol = (int) round((intersectPts[minIndex].getX() - leftX)/dx);
            }

            if(intersectPts[maxIndex].getX() >= rightX) {
                rightIntersectCol = screenWidth - 1;
            }
            else {
                rightIntersectCol = screenWidth - (1 + ((int) round((rightX - intersectPts[maxIndex].getX())/dx)));
            }

            // determining Za & Zb values
            Point p1 = tri_angle[i].corners[(int) intersectPts[minIndex].getZ()];
            Point p2 = tri_angle[i].corners[(int) intersectPts[minIndex].getW()];

            double tmp = (intersectPts[minIndex].getY() - p1.getY())*(p2.getZ() - p1.getZ())/(p2.getY() - p1.getY());
            double za = p1.getZ() + tmp;

            p1 = tri_angle[i].corners[(int) intersectPts[maxIndex].getZ()];
            p2 = tri_angle[i].corners[(int) intersectPts[maxIndex].getW()];

            double tmp2 = (intersectPts[maxIndex].getY() - p1.getY())*(p2.getZ() - p1.getZ())/(p2.getY() - p1.getY());
            double zb = p1.getZ() + tmp2;

            // scanning from leftIntersectCol to rightIntersectCol (incl)
            double zp;
            double consTerm = dx*(zb - za)/(intersectPts[maxIndex].getX() - intersectPts[minIndex].getX());
            double temp = ((leftX + leftIntersectCol*dx) - intersectPts[minIndex].getX())*(zb - za)/(intersectPts[maxIndex].getX()- intersectPts[minIndex].getX());

            for(int col = leftIntersectCol; col <= rightIntersectCol; col++) {
                // calculating z value
                if(col == leftIntersectCol) {
                    zp = za + temp;
                }
                else {
                    zp = zp + consTerm;
                }

                // comparing computed z value with current value in z_Buff[row][col] & frontLimitZ and updating accordingly
                if(zp > frontLimitZ && zp < z_Buff[row][col]) {
                    z_Buff[row][col] = zp;

                    frame_Buff[row][col].redVal = tri_angle[i].rgb.redVal;
                    frame_Buff[row][col].greenVal = tri_angle[i].rgb.greenVal;
                    frame_Buff[row][col].blueVal = tri_angle[i].rgb.blueVal;
                }
            }
        }
    }
}


void imageGeneration(){

    bitmap_image bitmapImg(screenWidth, screenHeight);

    for(int row = 0; row < screenHeight; row++) {
        for(int col = 0; col < screenWidth; col++) {
            bitmapImg.set_pixel(col, row, frame_Buff[row][col].redVal, frame_Buff[row][col].greenVal, frame_Buff[row][col].blueVal);
        }
    }

    bitmapImg.save_image("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/out.bmp");

    cout << "Image generated!" << endl;
}


void z_bufferWrite(ofstream & output){

    for(int row = 0; row < screenHeight; row++) {
        for(int col = 0; col < screenWidth; col++) {
            if(z_Buff[row][col] < rearLimitZ) {
                output << z_Buff[row][col] << '\t';
            }
        }
        output << endl;
    }
    output.close();

    cout << "z-buffer completed!" << endl;
}


void freeMemory(){

    for(int i=0; i<screenHeight; i++) {
        delete[] z_Buff[i];
    }
    delete[] z_Buff;

    for(int i=0; i<screenHeight; i++) {
        delete[] frame_Buff[i];
    }
    delete[] frame_Buff;

    cout << "Memory freed!" << endl;
}


int main(int argc, char** argv) {

    ifstream input("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/scene.txt",ios::in);
    ifstream input2("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/stage1.txt",ios::in);
    ifstream input3("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/stage2.txt",ios::in);
    ifstream input4("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/stage3.txt",ios::in);
    ifstream input5("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/config.txt",ios::in);

    ofstream output("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/stage1.txt",ios::out);
    ofstream output2("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/stage2.txt",ios::out);
    ofstream output3("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/stage3.txt",ios::out);
    ofstream output4("G:/Chrome Downloads/L-4 T-1 (2022)/CSE 410/Practice/Raster/test-cases/1/z-buffer.txt",ios::out);

    readScene(input);
    modelTransForm(input,output);

    viewTransForm(input2,output2);

    projectionTransForm(input3,output3);

    Triangle triangles[triagleCount];

    for(int i=0; i<triagleCount; i++)
        tri_angle.push_back(triangles[i]);

    readConfig(input5);
    readStage3(input4);
    initZ_buffer();
    applyProcedure();
    imageGeneration();
    z_bufferWrite(output4);
    freeMemory();

    cout << "Done!" << endl;
    return 0;
}
