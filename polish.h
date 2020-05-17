#ifndef __polish__
#define __polish__

#include <vector>
#include <iostream>
#include <cstring>

#define init_Layer 1
#define n_Layers 4
#define AR_min .7
#define AR_max 1.3

class Cnet;

class Cnets;

class Cmodule {


public:

    std::string m_name;
	std::vector<Cnet*> conn_net;
    char m_tag;
    int m_ptr;
    int m_parent;
    int m_left;
    int m_right;
    double m_x;
    double m_y;
    double m_z;
    double m_w;
    double m_d;
    double m_h;

/*The term t is the tag information, and its values are X, Y, Z for internal nodes and module names for leaves of the tree. The terms p, l, r denote cell locations of the parent, the left child and the right child of the
node. The terms x, y, z, w, d, h are the dimensional information of the corresponding module */

    Cmodule();
    Cmodule(int p, double w, double d, double h, std::string bname, int ptr );
    Cmodule(char t, int p, int left, int right);
    ~Cmodule();

    void print ();

    Cmodule& operator = (Cmodule&);
};


class Cmodules {

public:
    std::vector <Cmodule*> m_modules;

    void init(char* filename);

    double volume();

    int nebrmove(int b, int c);
    void perturb();


    void compute_vol (int module_id);
    double compute_wire(Cnets& my_cnets);
    int compute_MIV(Cnets& my_cnets);

    void coordinate(int module_id);
    void plot(const char* file_name);


    Cmodules();
    ~Cmodules();

    void print ();

    Cmodules& operator = (Cmodules&);



};

class Cnet{
public:
    int m_num;
    std::vector<Cmodule*>m_net;
    Cnet();
    ~Cnet();
    void print();
};

class Cnets{
public:
    std::vector<Cnet*>m_nets;

    Cnets();
    ~Cnets();
   void print();
};


#endif
