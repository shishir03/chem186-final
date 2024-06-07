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
    std::random_device rd;
    std::mt19937 gen(rd());

    System* s = new System(25, 300);

    std::normal_distribution<double> gaussian(0, 1);
    for(int i = 0; i < 25; i++) {
        // Initialize a new He atom at random coordinates
        double x = (rand() % 2500) / 100.0 - 12.5;
        double y = (rand() % 2500) / 100.0 - 12.5;
        double z = (rand() % 2500) / 100.0 - 12.5;

        double vx = sqrt(kB*300 / 2)*gaussian(gen);
        double vy = sqrt(kB*300 / 2)*gaussian(gen);
        double vz = sqrt(kB*300 / 2)*gaussian(gen);

        atom* a = new atom(x, y, z, vx, vy, vz, 0, 2);
        molecule* mol = new molecule();
        mol->atoms.push_back(a);
        s->molecules.push_back(mol);
    }

    s->run(10);
}