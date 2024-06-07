#include "molecule.h"
#include <tuple>

class System {
    public:
        double time = 0;
        std::vector<molecule*> molecules;
        // (epsilon, sigma)
        std::vector<std::pair<double, double>> lennard_jones_pots;

        System();
        double potential_energy();
        std::tuple<double, double, double> force(double init_pot_energy, atom* a);
        void do_timestep();
        void run(int num_steps);
};