#include "molecule.h"
#include <math.h>

atom::atom(double x, double y, double z, double charge, int mass) {
    this->x = x;
    this->y = y;
    this->z = z;
    
    vx = 0;
    vy = 0;
    vz = 0;

    this->charge = charge;
    this->mass = mass;
}

double atom::get_distance(atom* a2) {
    return sqrt((x - a2->x)*(x - a2->x) + (y - a2->y)*(y - a2->y) + (z - a2->z)*(z - a2->z));
}

double bond::get_length() {
    return a1->get_distance(a2);
}

bond::bond(atom* a1, atom* a2) {
    this->a1 = a1;
    this->a2 = a2;
}