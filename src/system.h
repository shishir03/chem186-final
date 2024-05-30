#include "molecule.h"
#include <tuple>

class System {
    public:
        double time = 0;
        std::vector<molecule*> molecules;

        double potential_energy(atom* a);
        std::tuple<double, double, double> force(atom* a);
        void do_timestep();
        void run(int num_steps);
};