#include "molecule.h"
#include <tuple>

class System {
    public:
        double time = 0;
        std::vector<molecule*> molecules;

        double potential_energy();
        std::tuple<double, double, double> force(double init_pot_energy, atom* a);
        void do_timestep();
        void run(int num_steps);
};