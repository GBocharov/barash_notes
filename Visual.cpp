#include <GL\glew.h>
#include <GL\freeglut.h>
#include "RigidBody.h"
#include "Visual.h"
#include <iostream>
#include <vector>


float last_time = 0.f;
RigidBody rigidBody = RigidBody();
RigidBody rigidBody_2 = RigidBody();
Vector g = { 0, G_constant, 0 };


RigidBody vec_to_RBody(std::vector<double> V)
{
    RigidBody b;
    b.coord = Vector{ V.at(0), V.at(1),V.at(2) };
    b.linear_moment = Vector{ V.at(3), V.at(4),V.at(5) };
    b.angular_moment = Vector{ V.at(6), V.at(7),V.at(8) };
    b.quaternion = Quaternion{ V.at(9), V.at(10), V.at(11),  V.at(12)};

    b.Inert_Tens_Inv = Matrix  { V.at(13), V.at(14),V.at(15), 
                        V.at(16), V.at(17),V.at(18), 
                        V.at(19), V.at(20),V.at(21),};

    b.Mass_Inv = V.at(22);

    return b;
}

std::vector<double> RBody_tovector(RigidBody b)
{
    std::vector<double> res;

    res.at(0) = b.coord.x;
    res.at(1) = b.coord.y;
    res.at(2) = b.coord.z;


    res.at(3) = b.linear_moment.x;
    res.at(4) = b.linear_moment.y;
    res.at(5) = b.linear_moment.z;

    res.at(6) = b.angular_moment.x;
    res.at(7) = b.angular_moment.y;
    res.at(8) = b.angular_moment.z;

    res.at(9) = b.quaternion.r;
    res.at(10) = b.quaternion.i;
    res.at(11) = b.quaternion.j;
    res.at(12) = b.quaternion.k;


    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            res.at(i + j) = b.Inert_Tens_Inv.values[i][j];

    res.at(22) = b.Mass_Inv;

    return res;
}

RigidBody f(RigidBody& body) { //D - func

    RigidBody result;

    result.coord = body.linear_moment / body.Mass_Inv; // dx

    result.Mass_Inv = body.Mass_Inv;

    body.quaternion = body.quaternion.normalize();

    body.Rotate_matrix = body.quaternion.toMatrix();

    Vector omega = ((body.Rotate_matrix * body.Inert_Tens_Inv) * body.Rotate_matrix.transpose()) * body.angular_moment;  // w

    result.quaternion = Quaternion{ 0.0, omega.x, omega.y, omega.z } * body.quaternion * 0.5;


    result.linear_moment = g * result.Mass_Inv; // t_F

    result.angular_moment = Vector{ 0, 0, 0 };      // t_torq

    result.Inert_Tens_Inv = body.Inert_Tens_Inv;

    return result;
}

template <typename T>
void solveRungeT(T& system, T(*f)(T& x), double h)
{
    T b1, b2, b3, b4;

    b1 = f(system);

    b2 = system + b1 * (h / 3);
    b2 = f(b2);

    b3 = system + ((b1 * (-h / 3)) + (b2 * h));
    b3 = f(b3);

    b4 = system + ((b1 * h) + (b2 * (-h)) + (b3 * h));
    b4 = f(b4);

    system = system + ((b1 * (1.0 / 8)) + (b2 * (3.0 / 8)) + (b3 * (3.0 / 8)) + (b4 * (1.0 / 8))) * h;
}

void solveBody(RigidBody& b, double h)
{
    solveRungeT(b, f, h);
    b.quaternion = b.quaternion.normalize();
    b.Rotate_matrix = b.quaternion.toMatrix();
}

void solveRunge(RigidBody& body, RigidBody(*f)(RigidBody& body), double h) {     
    
    RigidBody b1(rigidBody.Inert_Tens_Inv), b2(rigidBody.Inert_Tens_Inv), b3(rigidBody.Inert_Tens_Inv), b4(rigidBody.Inert_Tens_Inv); // I_Inv = const 

    b1 = f(body);

    b2 = body + b1 * (h / 3);
    b2 = f(b2);

    b3 = body + ((b1 * (-h / 3)) + (b2 * h));
    b3 = f(b3);

    b4 = body + ((b1 * h) + (b2 * (-h)) + (b3 * h));
    b4 = f(b4);

    body = body + ((b1 * (1.0 / 8)) + (b2 * (3.0 / 8)) + (b3 * (3.0 / 8)) + (b4 * (1.0 / 8))) * h;

    body.quaternion = body.quaternion.normalize();

    body.Rotate_matrix = body.quaternion.toMatrix();
}




void solve(RigidBody& body, double step, double e)
{
    double old_step = 0;
    double currentStep = step;
    Vector collisionPlacement = {0,0,0};
    collisionPlacement = { 0, 0, 0};
    RigidBody next_body = body;

    solveBody(next_body, step);
    RigidBody next_next_body = next_body;

    solveBody(next_next_body, step);

    if (body.check_collision() && (! next_body.check_collision()  ))
    {
       // body.print_Spheres();
        solveBody(body, step);
        std::cout << "----------------------------------->kostyl" << std::endl;
        return;
    }


    

    if (next_body.check_collision() && next_next_body.check_collision())
    {
        int k = 0;
        while (1)
        {
            

            next_body = body;
            solveBody(next_body, currentStep);

            //next_body.print_Spheres();
            bool is_collision_defined = next_body.find_collision_place(e, collisionPlacement);

            if (is_collision_defined) 
            {
                break;
            }
            else
            {
                if (next_body.check_collision())
                {
                    double t = currentStep;
                    currentStep = currentStep - abs(old_step - currentStep) / 2;
                    old_step = t;
                }
                else 
                {
                    double t = currentStep;
                    currentStep = currentStep + abs(old_step - currentStep) / 2;
                    old_step = t;
                }    
            }

        }
        solveBody(body, currentStep);
        body.floorCollisionAction(collisionPlacement);
        solve(body, step - currentStep, e);
        
    }
    else 
    {
        solveBody(body, step);
    }
}


void Idle() {
    glutPostRedisplay();
}

void drawFloor()
{
    glPushMatrix();
    glTranslated(0,-20,0);

    glBegin(GL_QUADS);
    glColor3d(0.5, 0.5, 0.5);
    glVertex3f(-200, FLOORLEVEL, -850);
    glVertex3f(200, FLOORLEVEL, -850);
    glVertex3f(200, FLOORLEVEL, -55);
    glVertex3f(-200, FLOORLEVEL, -55);

    glEnd();

    glPopMatrix();
}

void drawS(Vector x, double r)
{
    glColor3f(1, 0, 0);

    glPushMatrix();

    glTranslated(x.x, x.y, x.z);

    glutWireSphere(r, 10, 10);

    glPopMatrix();
}

void drawLINE(Vector x, Vector y)
{
    glColor3f(0, 1, 0);

   // glPushMatrix();

    glLineWidth(3); 
    glBegin(GL_LINES); 
    glVertex3d(x.x, x.y, x.z); 
    glVertex3d(y.x, y.y, y.z);
    glEnd();

   // glPopMatrix();
}

void drawBody() {

    for (int i = 0; i < rigidBody.Spheres.size(); i++) {
        drawS(rigidBody.Spheres.at(i).coord, rigidBody.Spheres.at(i).Radius);

        if (i < rigidBody.Spheres.size() - 1)
            drawLINE(rigidBody.Spheres.at(i).coord, rigidBody.Spheres.at(i + 1).coord);
        else
            drawLINE(rigidBody.Spheres.at(i).coord, rigidBody.Spheres.at(1).coord);
    }

    //glColor3f(0, 1, 0);
    //glutWireSphere(50, 20, 20);
}

void drawBody1() {

    for (int i = 0; i < rigidBody.Spheres.size(); i++) {
        drawS(rigidBody.Spheres.at(i).coord, rigidBody.Spheres.at(i).Radius);

        if (i < rigidBody.Spheres.size() - 1)
            drawLINE(rigidBody.Spheres.at(i).coord, rigidBody.Spheres.at(i + 1).coord);
    }

    //glColor3f(0, 1, 0);
    //glutWireSphere(50, 20, 20);
}

void initBody() 
{
    rigidBody.initSpheres1();
    rigidBody.findCenter();
    rigidBody.findTensor();
    rigidBody.print();
}


float log_time = 0;

void Display() {
    glViewport(0, 0, 600, 600);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 1, 500);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    float now = clock() / (float)CLOCKS_PER_SEC;
    float dt = now - last_time;
    last_time = now;
    log_time += dt;

    solve(rigidBody, dt, 0.1);   //  solve (body, dt,  epsilon) //0.000005

    if (log_time >= 1) {
        //std::cout << "E = " << rigidBody.getEnergyInv() << std::endl;
        log_time = 0;
    }
    glPushMatrix();

    drawFloor();

    glTranslated(rigidBody.coord.x, rigidBody.coord.y, rigidBody.coord.z - 300);

    glRotated(acos(rigidBody.quaternion.r) * 360 / 3.14, rigidBody.quaternion.i, rigidBody.quaternion.j, rigidBody.quaternion.k);

    drawBody1();

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}


int main(int argc, char* argv[]) {
    glutInit(&argc, argv);

    initBody();

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("W_M");
    glutDisplayFunc(Display);
    glutIdleFunc(Idle);
    glEnable(GL_DEPTH_TEST);
    glutMainLoop();
    return 0;
}