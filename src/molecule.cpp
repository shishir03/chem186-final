#include "molecule.h"
#include <math.h>

double bond::get_length() {
    return sqrt((a1->x - a2->x)*(a1->x - a2->x) + (a1->y - a2->y)*(a1->y - a2->y) + (a1->z - a2->z)*(a1->z - a2->z));
}

bond::bond(atom* a1, atom* a2) {
    this->a1 = a1;
    this->a2 = a2;
}