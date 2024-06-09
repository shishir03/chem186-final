#include <vector>

class atom {
    public:
        double x;
        double y;
        double z;
        double charge;
        int mass;
        int atomic_no;

        double vx;
        double vy;
        double vz;

        atom(double x, double y, double z, double vx, double vy, double vz, double charge, int atomic_no, int mass);
        double get_distance(atom* a2);
        double get_speed();
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
