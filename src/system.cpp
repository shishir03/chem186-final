#include "system.h"
#include "consts.h"
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <sstream>

// Initialize Lennard Jones potentials array
System::System() {
    std::ifstream file("lennard_jones.txt");

    std::string line;

    // Read the file line by line
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string value;

        // Parse each line by commas
        while (std::getline(ss, value, ',')) {
            row.push_back(value);
        }

        std::pair<double, double> vals = std::make_pair(std::stod(row[3]), std::stod(row[4]));

        // Add the row to the data
        lennard_jones_pots.push_back(vals);
    }

    // Close the file
    file.close();
}

// Units of eV A^-1
std::tuple<double, double, double> System::force(double init_pot_energy, atom* a) {
    a->x += dr;
    double new_pot_energy = potential_energy();
    double fx = -(new_pot_energy - init_pot_energy) / dr;
    a->x -= dr;

    a->y += dr;
    new_pot_energy = potential_energy();
    double fy = -(new_pot_energy - init_pot_energy) / dr;
    a->y -= dr;

    a->z += dr;
    new_pot_energy = potential_energy();
    double fz = -(new_pot_energy - init_pot_energy) / dr;
    a->z -= dr;

    return std::make_tuple(fx, fy, fz);
}

// Units of eV
double System::potential_energy() {
    double total_potential = 0;

    for(auto mol : molecules) {
        for(auto b : mol->bonds) {
            // Bond stretch term
            double l = b->get_length();
            double dl = l - equil;
            total_potential += 0.5*k*dl*dl;
        }

        for(auto a : mol->atoms) {
            for(auto mol2 : molecules) {
                for(auto a2 : mol2->atoms) {
                    if(a != a2) {
                        double r = a->get_distance(a2);
                        // Lennard Jones potential
                        std::pair<double, double> pot1 = lennard_jones_pots[a->mass - 1];
                        std::pair<double, double> pot2 = lennard_jones_pots[a2->mass - 1];
                        // units of eV
                        double total_epsilon = sqrt(std::get<0>(pot1)*std::get<0>(pot2));
                        double total_sigma = (std::get<1>(pot1) + std::get<1>(pot2)) / 2;

                        total_potential += 4*total_epsilon*(pow(total_sigma / r, 12) - pow(total_sigma / r, 6));
                    }
                }
            }
        }
    }

    return total_potential;
}

void System::do_timestep() {
    for(auto mol : molecules) {
        for(auto a : mol->atoms) {
            double u = potential_energy();

            // Compute forces
            std::tuple<double, double, double> Fi = force(u, a);
            // Acceleration: units of eV A^-1 amu^-1 => 9648.53322 A ps^-2
            double ax = accel_constant * std::get<0>(Fi) / a->mass;
            double ay = accel_constant * std::get<1>(Fi) / a->mass;
            double az = accel_constant * std::get<2>(Fi) / a->mass;

            // printf("Acceleration: %.2f %.2f %.2f\n", ax, ay, az);

            // Update positions
            a->x += a->vx*dt + 0.5*ax*dt*dt;
            a->y += a->vy*dt + 0.5*ay*dt*dt;
            a->z += a->vz*dt + 0.5*az*dt*dt;

            u = potential_energy();

            // Compute new forces
            std::tuple<double, double, double> Fi1 = force(u, a);
            double ax1 = std::get<0>(Fi1) / a->mass;
            double ay1 = std::get<1>(Fi1) / a->mass;
            double az1 = std::get<2>(Fi1) / a->mass; 

            // Update velocities
            a->vx += 0.5*(ax + ax1)*dt;
            a->vy += 0.5*(ay + ay1)*dt;
            a->vz += 0.5*(az + az1)*dt;
        }
    }

    time += dt;
}

void System::run(int num_timesteps) {
    for(int i = 0; i < num_timesteps; i++) {
        printf("Timestep %d\n", i);
        for(int j = 0; j < molecules.size(); j++) {
            printf("Molecule %d ", j);
            for(auto atom : molecules[j]->atoms) {
                printf("(%.5f, %.5f, %.5f)\n", atom->x, atom->y, atom->z);
            }
        }

        do_timestep();
    }
}
