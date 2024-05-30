#include "system.h"
#include "consts.h"
#include <math.h>
#include <stdio.h>

std::tuple<double, double, double> System::force(atom* a) {
    std::vector<bond*> bonds = ((molecule*) a->mol)->bonds;
    for(auto b : bonds) {
        if(b->a1 == a) {
            double l = b->get_length();
            double dl = l - eq;

            double dx = dl*(b->a1->x - b->a2->x) / l;
            double dy = dl*(b->a1->y - b->a2->y) / l;
            double dz = dl*(b->a1->z - b->a2->z) / l;

            return std::make_tuple(-k*dx, -k*dy, -k*dz);
        } else if(b->a2 == a) {
            double l = b->get_length();
            double dl = l - eq;

            double dx = dl*(b->a2->x - b->a1->x) / l;
            double dy = dl*(b->a2->y - b->a1->y) / l;
            double dz = dl*(b->a2->z - b->a1->z) / l;

            return std::make_tuple(-k*dx, -k*dy, -k*dz);
        }
    }

    return std::make_tuple(0, 0, 0);
}

double System::potential_energy(atom* a) {
    return 0;
}

void System::do_timestep() {
    for(auto mol : molecules) {
        for(auto a : mol->atoms) {
            // Compute forces
            std::tuple<double, double, double> Fi = force(a);
            double ax = std::get<0>(Fi) / a->mass;
            double ay = std::get<1>(Fi) / a->mass;
            double az = std::get<2>(Fi) / a->mass;

            printf("Acceleration: %.2f %.2f %.2f\n", ax, ay, az);

            // Update positions
            a->x += a->vx*dt + 0.5*ax*dt*dt;
            a->y += a->vy*dt + 0.5*ay*dt*dt;
            a->z += a->vz*dt + 0.5*az*dt*dt;

            // Compute new forces
            std::tuple<double, double, double> Fi1 = force(a);
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
