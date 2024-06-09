#include "molecule.h"
#include "consts.h"
#include <tuple>
#include <random>

class System {
    public:
        double t = 0;
        std::vector<molecule*> molecules;
        // (epsilon, sigma)
        std::vector<std::pair<double, double>> lennard_jones_pots;
        double box_size;
        double temper;

        System(double size, double temper);
        void sample_boltzmann(atom* a);
        void run(int num_steps);

    private:
        double kinetic_energy();
        double potential_energy();
        double temperature();
        void run_thermostat();
        std::tuple<double, double, double> force(double init_pot_energy, atom* a);
        void do_timestep();
};