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
#include <random>

int main() {
    srand(time(NULL));
    std::random_device rd;
    std::mt19937 gen(rd());

    double box_size = 10;
    double temper = 300;
    double atomic_no = 2;   // Use He atoms for now
    double atomic_mass = 4;

    System* s = new System(box_size, temper);

    std::normal_distribution<double> gaussian(0, 1);
    for(int i = 0; i < 25; i++) {
        // Initialize a new atom at random coordinates
        double x = (rand() % ((int)box_size * 100)) / 100.0 - box_size/2;
        double y = (rand() % ((int)box_size * 100)) / 100.0 - box_size/2;
        double z = (rand() % ((int)box_size * 100)) / 100.0 - box_size/2;

        double vx = sqrt(kB*300 / atomic_mass)*gaussian(gen);
        double vy = sqrt(kB*300 / atomic_mass)*gaussian(gen);
        double vz = sqrt(kB*300 / atomic_mass)*gaussian(gen);

        atom* a = new atom(x, y, z, vx, vy, vz, 0, atomic_no, atomic_mass);
        molecule* mol = new molecule();
        mol->atoms.push_back(a);
        s->molecules.push_back(mol);
    }

    s->run(10);
}