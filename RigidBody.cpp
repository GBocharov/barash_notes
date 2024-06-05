#include "RigidBody.h"
#include <iostream>


RigidBody RigidBody::operator+(RigidBody A) const {

    RigidBody result;
    result.coord = coord + A.coord;
    result.quaternion = quaternion + A.quaternion;
    result.linear_moment = linear_moment + A.linear_moment;
    result.angular_moment = angular_moment + A.angular_moment;

    result.Spheres = Spheres;

    result.Inert_Tens_Inv = Inert_Tens_Inv;

    result.Mass_Inv = Mass_Inv;

    return result;
}

RigidBody RigidBody::operator*(double h) const {

    RigidBody result{};

    result.coord = coord * h;

    result.quaternion = quaternion * h;

    result.linear_moment = linear_moment * h;

    result.angular_moment = angular_moment * h;

    result.Spheres = Spheres;

    result.Inert_Tens_Inv = Inert_Tens_Inv;

    result.Mass_Inv = Mass_Inv;

    return result;
}

RigidBody::RigidBody() {
  
    coord = { 0,0,0 };
    quaternion = { cos(0), 0, 0, 0 };

    Rotate_matrix = quaternion.normalize().toMatrix();
   
    linear_moment = { 0,0,0 };
    angular_moment = Vector{ 10, 10, 10 }; // { {0, 10000, 0} }

    //I_Inv.printMatr();

}

RigidBody::RigidBody(Matrix m) {
    Inert_Tens_Inv = m;
    //I_Inv.printMatr();
}


void RigidBody::initSpheres1()
{
    double m = 0.1;
    double dist_k = 10;
    double fix_radius = 7;

    Vector x1 = Vector{ 2, 2, 2 }*dist_k;
    double R1 = fix_radius;
    double mass1 = m;
    Sphere s1 = Sphere(x1, mass1, R1);


    Vector x2 = Vector{ 0, 0, 0 }*dist_k ;
    double R2 = fix_radius;
    double mass2 = m;
    Sphere s2 = Sphere(x2, mass2, R2);

    Vector x3 = Vector{ 0, 0, 4 }*dist_k;
    double R3 = fix_radius;
    double mass3 = m;
    Sphere s3 = Sphere(x3, mass3, R3);

    Vector x4 = Vector{ 4, 0, 0 }*dist_k;
    double R4 = fix_radius;
    double mass4 = m;
    Sphere s4 = Sphere(x4, mass4, R4);


    Vector x5 = Vector{ 4, 0, 4 }*dist_k;
    double R5 = fix_radius;
    double mass5 = m;
    Sphere s5 = Sphere(x5, mass5, R5);


    Vector x6 = Vector{ 0, 4, 0 }*dist_k;
    double R6 = fix_radius;
    double mass6 = m;
    Sphere s6 = Sphere(x6, mass6, R6);

    Vector x7 = Vector{ 0, 4, 4 }*dist_k;
    double R7 = fix_radius;
    double mass7 = m;
    Sphere s7 = Sphere(x7, mass7, R7);

    Vector x8 = Vector{ 4, 4, 0 }*dist_k;
    double R8 = fix_radius;
    double mass8 = m;
    Sphere s8 = Sphere(x8, mass8, R8);


    Vector x9 = Vector{ 4, 4, 4 }*dist_k;
    double R9 = fix_radius;
    double mass9 = m;
    Sphere s9 = Sphere(x9, mass9, R9);

    // test1 cube
   {

      Spheres.push_back(s1);
      Spheres.push_back(s2);
      Spheres.push_back(s3);
      Spheres.push_back(s4);
      Spheres.push_back(s5);
      Spheres.push_back(s6);
      Spheres.push_back(s7);
      Spheres.push_back(s8);
      Spheres.push_back(s9);/**/

    }
   /*
    {
        Spheres.push_back(s2);
        Spheres.push_back(s3);
        Spheres.push_back(s5);
        Spheres.push_back(s6);
        Spheres.push_back(s7);
        Spheres.push_back(s9);
    }*/
}

void RigidBody::findCenter() 
{
    double total_mass = 0;
    Vector x = {0, 0, 0};

    for (auto &i : Spheres) 
    {
        x = x + i.coord * i.mass;
        total_mass += i.mass;
    }

    x = x / total_mass;

    Mass_Inv = total_mass;

    std::cout << "Center mass = ";
    x.print();

    coord = x;

    for (auto& i : Spheres) 
    {
        i.coord = i.coord - x;
    }

}


void RigidBody::findTensor() 
{

    for (auto& i : Spheres) 
    {
       // std::cout << "DO:" << std::endl;
       // i.I.printMatr();

        i.Compute({ 0.0, 0.0, 0.0 });

      //  std::cout << "Posle:" << std::endl;
      //  i.I.printMatr();
    }

    for (auto& i : Spheres) 
    {
        Inert_Tens_Inv = Inert_Tens_Inv + i.I;
    }

    std::cout << "Final inertia tensor:" << std::endl;

    Inert_Tens_Inv.printMatr();

    Inert_Tens_Inv = Inert_Tens_Inv.Inverse();
    
    std::cout << "Final inertia tensor inverse:" << std::endl;

    Inert_Tens_Inv.printMatr();
}


double RigidBody::getEnergyInv() 
{
    double g = abs(G_constant);                 // |g| !
    double h = coord.y - FLOORLEVEL;
    double E_p = Mass_Inv * g * h;

    double v = linear_moment.abs() / Mass_Inv;
    double E_linear = (Mass_Inv * v * v) / 2;

    Vector w = (Rotate_matrix * Inert_Tens_Inv * Rotate_matrix.transpose()) * angular_moment;

    double E_angular = (w * (Rotate_matrix * Inert_Tens_Inv * Rotate_matrix.transpose()).Inverse()).dot(w) / 2;
    std::cout << "E_p " << E_p << "     E_linear " << E_linear << "     E_angular " << E_angular << std::endl;
    return (E_p + E_linear + E_angular);
}

void RigidBody::print_Spheres() 
{
    for (auto& i : Spheres) 
    {
        std::cout << "Center : ";
        (Rotate_matrix * i.coord + coord).print();

    }
}


void RigidBody::print()
{
    std::cout << "-------------------" << std::endl;
    std::cout << "Mass = "<< Mass_Inv<<std::endl;
    std::cout << "x = ";
    coord.print();
    std::cout << "l = ";

    linear_moment.print();
    std::cout << "L = ";

    angular_moment.print();


    std::cout << "I_Inv = "<<std::endl;

    Inert_Tens_Inv.printMatr();


    std::cout << "-------------------" << std::endl;
}


void RigidBody::floorCollisionAction(Vector collisionPlace)
{
    Vector coll_Pl_Rotate = Rotate_matrix * collisionPlace; 
    Matrix inertiaTensorInv = Rotate_matrix * Inert_Tens_Inv * Rotate_matrix.transpose();
    Vector p_dot = linear_moment / Mass_Inv + ((inertiaTensorInv * angular_moment)* coll_Pl_Rotate);

    Vector n;
    n = { 0, 1, 0 };
    double numerator = (-(1 + 1) * n.dot(p_dot));
    double term1 = 1 / Mass_Inv;
    double term2 = n.dot((inertiaTensorInv * (coll_Pl_Rotate *(n)))*(coll_Pl_Rotate));
    Vector total_force =   n*( (numerator / (term1 + term2)));

    linear_moment = linear_moment + total_force;
    angular_moment = angular_moment + coll_Pl_Rotate * (total_force);

    std::cout << "-------<COLLISION ACTION START>-----" << std::endl;

    std::cout << "      <Energy>     " << std::endl;;
    std::cout << std::endl << "total energy = " << getEnergyInv() << std::endl;
    std::cout << "    <COLLISION PLACE>    "<< std::endl;
    (coll_Pl_Rotate + coord + Vector{0,-7,0}).print();

    std::cout << "-------<COLLISION ACTION END>-----"<<std::endl;
}






//--------------------------------------------------NEW-------------------------------------------------------------

bool RigidBody::check_collision ()
{
    double collideTolerance = COLLIDETOLERANCE;
    double floorLevel = FLOORLEVEL;
    Vector floor_normal = { 0, 1, 0 };   // floar normal

   Rotate_matrix = quaternion.normalize().toMatrix();

    for (auto& i : Spheres)
    {
        Vector curent_coord = Rotate_matrix * i.coord;
        Matrix inertiaTensorInv = Rotate_matrix * Inert_Tens_Inv * Rotate_matrix.transpose();
        Vector p_dot = linear_moment / Mass_Inv + ((inertiaTensorInv * angular_moment) * (curent_coord));

        bool isBelowFloor = ((curent_coord + coord).dot(floor_normal) - (floorLevel + i.Radius )  < 0.01);

        bool isGoesToFloor = p_dot.dot(floor_normal) <= -0.001; // -collideTolerance;

        if (isBelowFloor && isGoesToFloor)
        {
            return true;
        }
    }

    return false;

}





bool RigidBody::find_collision_place(double e, Vector& res)
{
    double floorLevel = FLOORLEVEL;
    double final_dist = 10000000000;
    int k = 0;
    double final_Radius = 0;
    Vector final_coord;
    int number = 0;
    int max_number = 0;



    for (auto& i : Spheres) // chose most collised sphere
    {
        number++;
        double distance = (Rotate_matrix * i.coord + coord).y - (floorLevel + i.Radius) ;

        Vector curent_coord = Rotate_matrix * i.coord;
        Vector floor_normal = { 0, 1, 0 };   // floar normal
        Vector p_dot = linear_moment / Mass_Inv + (((Rotate_matrix * Inert_Tens_Inv * Rotate_matrix.transpose()) * angular_moment) * (curent_coord));
        bool isGoesToFloor = p_dot.dot(floor_normal) <= 0; // -collideTolerance;

        if (distance < final_dist && isGoesToFloor)
        {
            final_dist = distance;
            final_coord = i.coord;
            final_Radius = i.Radius;
            max_number = number;

        }

    }

    if (abs(final_dist) < e )
    {
        k = 1; // find point (Sphere center)
    }
    if (k == 1) {
        // std::cout << "Collision coord_W --> " << std::endl; // collision coord 
         //((Rotate_matrix * final_coord + coord + Vector{ 0, -final_Radius, 0 })).print(); // vector in W.coord
        //res = (final_coord + Rotate_matrix.transpose() * Vector { 0, -final_Radius, 0 });
        res = (final_coord);
        //std::cout << "dist == " << final_dist << std::endl;
        return true;
    }
    return false;
}




void RigidBody::initSpheres2()
{
    double m = 0.1;
    double dist_k = 10;
    double fix_radius = 7;

    Vector x1 = Vector{ 0, 0, 0 }*dist_k;
    double R1 = fix_radius;
    double mass1 = m;
    Sphere s1 = Sphere(x1, mass1, R1);


    Vector x2 = Vector{ 1, 0, 0 }*dist_k;
    double R2 = fix_radius;
    double mass2 = m;
    Sphere s2 = Sphere(x2, mass2, R2);

    Vector x3 = Vector{ 2, 0, 4 }*dist_k;
    double R3 = fix_radius;
    double mass3 = m;
    Sphere s3 = Sphere(x3, mass3, R3);

    Vector x4 = Vector{ 3, 0, 0 }*dist_k;
    double R4 = fix_radius;
    double mass4 = m;
    Sphere s4 = Sphere(x4, mass4, R4);


    Vector x5 = Vector{ 2, -1, 0 }*dist_k;
    double R5 = fix_radius;
    double mass5 = m;
    Sphere s5 = Sphere(x5, mass5, R5);


    Vector x6 = Vector{ 1, -2, 0 }*dist_k;
    double R6 = fix_radius;
    double mass6 = m;
    Sphere s6 = Sphere(x6, mass6, R6);

    Vector x7 = Vector{ 0, -3, 0 }*dist_k;
    double R7 = fix_radius;
    double mass7 = m;
    Sphere s7 = Sphere(x7, mass7, R7);

    Vector x8 = Vector{ 1, -3, 0 }*dist_k;
    double R8 = fix_radius;
    double mass8 = m;
    Sphere s8 = Sphere(x5, mass5, R5);


    Vector x9 = Vector{ 2, -3, 0 }*dist_k;
    double R9 = fix_radius;
    double mass9 = m;
    Sphere s9 = Sphere(x6, mass6, R6);

    Vector x10 = Vector{ 3, -3, 0 }*dist_k;
    double R10 = fix_radius;
    double mass10 = m;
    Sphere s10 = Sphere(x7, mass7, R7);


    // test1 cube
    {

        Spheres.push_back(s1);
        Spheres.push_back(s2);
        Spheres.push_back(s3);
        Spheres.push_back(s4);
        Spheres.push_back(s5);
        Spheres.push_back(s6);
        Spheres.push_back(s7);
        Spheres.push_back(s8);
        Spheres.push_back(s9);
        Spheres.push_back(s10);/**/

    }
    /*
     {
         Spheres.push_back(s2);
         Spheres.push_back(s3);
         Spheres.push_back(s5);
         Spheres.push_back(s6);
         Spheres.push_back(s7);
         Spheres.push_back(s9);
     }*/
}