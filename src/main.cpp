/**
 * Plan:
 *  - get basic velocity verlet method working
 *  - figure out how to get AMBER forcefield working (either via a library or doing it manually)
 *  - add support for parsing mol2 files to get molecule data
 *  - figure out how to display molecules (if there's time)
*/
#include "system.h"

int main() {
    System* s = new System();
    atom a1 = { 0, 0, 0, 0, 2, 0, 0, 0, nullptr };
    atom a2 = { 5, 0, 0, 0, 2, 0, 0, 0, nullptr };
    bond* b = new bond(&a1, &a2);

    molecule* mol = new molecule();
    mol->bonds.push_back(b);
    mol->atoms.push_back(&a1);
    mol->atoms.push_back(&a2);

    a1.mol = mol;
    a2.mol = mol;

    s->molecules.push_back(mol);

    s->run(50);
}