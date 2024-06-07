#include "molecule.h"
#include "consts.h"
#include <tuple>

class System {
    public:
        double time = 0;
        std::vector<molecule*> molecules;
        // (epsilon, sigma)
        std::vector<std::pair<double, double>> lennard_jones_pots;
        double box_size;
        double temper;

        System(double size, double temper);
        void run(int num_steps);

    private:
        double potential_energy();
        std::tuple<double, double, double> force(double init_pot_energy, atom* a);
        void do_timestep();
};