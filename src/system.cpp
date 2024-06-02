#include "system.h"
#include "consts.h"
#include <math.h>
#include <stdio.h>

// This can probably be parallelized
std::tuple<double, double, double> System::force(double init_pot_energy, atom* a) {
    a->x += dr;
    double new_pot_energy = potential_energy();
    double fx = -(new_pot_energy - init_pot_energy) / dr;
    a->x -= dr;
    // printf("%.2f %.2f\n", init_pot_energy, new_pot_energy);

    a->y += dr;
    new_pot_energy = potential_energy();
    double fy = -(new_pot_energy - init_pot_energy) / dr;
    a->y -= dr;
    // printf("%.2f %.2f\n", init_pot_energy, new_pot_energy);

    a->z += dr;
    new_pot_energy = potential_energy();
    double fz = -(new_pot_energy - init_pot_energy) / dr;
    a->z -= dr;
    // printf("%.2f %.2f\n", init_pot_energy, new_pot_energy);

    return std::make_tuple(fx, fy, fz);
}

double System::potential_energy() {
    double total_potential = 0;

    for(auto mol : molecules) {
        for(auto b : mol->bonds) {
            double l = b->get_length();
            double dl = l - eq;
            total_potential += 0.5*k*dl*dl;
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
            double ax = std::get<0>(Fi) / a->mass;
            double ay = std::get<1>(Fi) / a->mass;
            double az = std::get<2>(Fi) / a->mass;

            printf("Acceleration: %.2f %.2f %.2f\n", ax, ay, az);

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
            printf("Molecule %d\n", j);
            for(auto atom : molecules[j]->atoms) {
                printf("Atom at %.2f %.2f %.2f\n", atom->x, atom->y, atom->z);
            }
        }

        do_timestep();
    }
}
