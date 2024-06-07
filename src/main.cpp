/**
 * Plan:
 *  - implement electrostatic forces
 *  - find library for bond stretch / equilibrium constants
 * If there's time:
 *  - figure out how to do bond torsion angles / dihedral angles
 *  - add support for parsing mol2 files to get molecule data
 *  - figure out how to display molecules
*/
#include "system.h"
#include <random>

int main() {
    srand(time(NULL));
    System* s = new System();

    for(int i = 0; i < 25; i++) {
        // Initialize a new He atom at random coordinates
        double x = (rand() % 200) / 100.0;
        double y = (rand() % 200) / 100.0;
        double z = (rand() % 200) / 100.0;

        atom* a = new atom(x, y, z, 0, 2);
        molecule* mol = new molecule();
        mol->atoms.push_back(a);
        s->molecules.push_back(mol);
    }

    s->run(10);
}