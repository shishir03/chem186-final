/**
 * Plan:
 *  - implement electrostatic forces
 *  - simple thermostat
 * If there's time:
 *  - find library for bond stretch / equilibrium constants
 *  - figure out how to do bond torsion angles / dihedral angles
 *  - add support for parsing mol2 files to get molecule data
 *  - figure out how to display molecules
*/
#include "system.h"

int main() {
    srand(time(NULL));

    double box_size = 10;
    double temper = 300;
    double atomic_no = 2;   // Use He atoms for now
    double atomic_mass = 4;

    System* s = new System(box_size, temper);

    for(int i = 0; i < 25; i++) {
        // Initialize a new atom at random coordinates
        double x = (rand() % ((int)box_size * 100)) / 100.0 - box_size/2;
        double y = (rand() % ((int)box_size * 100)) / 100.0 - box_size/2;
        double z = (rand() % ((int)box_size * 100)) / 100.0 - box_size/2;

        atom* a = new atom(x, y, z, 0, 0, 0, 0, atomic_no, atomic_mass);
        s->sample_boltzmann(a);

        molecule* mol = new molecule();
        mol->atoms.push_back(a);
        s->molecules.push_back(mol);
    }

    s->run(10);
}