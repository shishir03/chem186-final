#include "molecule.h"
#include <math.h>

atom::atom(double x, double y, double z, double vx, double vy, double vz, double charge, int atomic_no, int mass) {
    this->x = x;
    this->y = y;
    this->z = z;

    this->vx = vx;
    this->vy = vy;
    this->vz = vz;

    this->charge = charge;
    this->atomic_no = atomic_no;
    this->mass = mass;
}

double atom::get_distance(atom* a2) {
    return sqrt((x - a2->x)*(x - a2->x) + (y - a2->y)*(y - a2->y) + (z - a2->z)*(z - a2->z));
}

double atom::get_speed() {
    return sqrt(vx*vx + vy*vy + vz*vz);
}

double bond::get_length() {
    return a1->get_distance(a2);
}

bond::bond(atom* a1, atom* a2) {
    this->a1 = a1;
    this->a2 = a2;
}