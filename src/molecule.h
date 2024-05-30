#include <vector>

struct atom {
    double x;
    double y;
    double z;
    double charge;
    int mass;

    double vx;
    double vy;
    double vz;

    void* mol;
};

class bond {
    public:
        atom* a1;
        atom* a2;
        double length;
        double angle;

        bond(atom* a1, atom* a2);
        double get_length();
};

struct molecule {
    std::vector<atom*> atoms;
    std::vector<bond*> bonds;
};
