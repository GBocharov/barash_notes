#pragma once

#include "Matrix.h"
#include "Quaternion.h"
#include "Vector.h"
#include "Sphere.h"
#include <vector>


#define G_constant -10
#define FLOORLEVEL -40

const double COLLIDETOLERANCE = -1.0;


 class RigidBody {
    public:
    std::vector<Sphere> Spheres;

    double Mass_Inv; // mass_invar

    Matrix Rotate_matrix{}; // rotate matrix
    Vector coord{}, linear_moment{}, angular_moment{}; // coord, l_momen, a_moment
    Quaternion quaternion{}; // quat
    Matrix Inert_Tens_Inv;   // inetr_tens_inversed

    RigidBody();
    RigidBody(Matrix m);
    RigidBody f();
    RigidBody operator+(RigidBody A) const;
    RigidBody operator*(double h) const;

    void initSpheres2();   
    void initSpheres1();  
    void findTensor();   
    void findCenter();  // find mass center
    void print();

    bool isCollided();

    bool findCollisionPlace(double e, Vector& res);
    void floorCollisionAction(Vector collisionPlace);
    bool check_collision();
    bool find_collision_place(double e, Vector& res);
    double getEnergyInv();
    void print_Spheres();

};
