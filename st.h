#ifndef __VSP__
#define __VSP__

#include <vector>
#include <iostream>
#include <cstring>

#define init_Layer 1
#define n_Layers 4
#define VSP_POS 0
#define VSP_NEG 1

#define VSP_NODE_VISITED 0
#define VSP_NODE_UNVISITED 1
#define VSP_NODE_TO_VISIT 2

class VSPGraphNode;
class VSPGraphEdge;  // directed edge
class VSPNet;


class VSPBlock {
  public:
    int m_id;
    std::string m_name;
    std::vector<VSPNet*> conn_net;

    double m_width;
    double m_depth;
    double m_height;
    double m_x;
    double m_y;
    double m_z;

    bool m_rotated;

  public:
    VSPBlock (int id);
    ~VSPBlock ();
};

class VSPNet{
public:
    int m_num;
    std::vector<VSPBlock*>m_net;
    VSPNet();
    ~VSPNet();

};

class VSPGraphNode {
  public:
    int m_id;

    std::vector<VSPGraphEdge*> m_edges_parents;  // incoming edges
    std::vector<VSPGraphEdge*> m_edges_children;  // outgoing edges

    int m_status;
    void* m_data1;
    double m_loc;

  public:
    VSPGraphNode (int id);
    ~VSPGraphNode ();
};



class VSPGraphEdge {
  public:
    VSPGraphNode* m_from;
    VSPGraphNode* m_to;

  public:
    VSPGraphEdge (VSPGraphNode* from, VSPGraphNode* to);
    ~VSPGraphEdge ();
};



class VSPGraph {
  public:
    std::vector<VSPGraphNode*> m_nodes;
    std::vector<VSPGraphEdge*> m_edges;

    VSPGraphNode* m_src_node;
    VSPGraphNode* m_tail_node;

  public:
    VSPGraph ();
    ~VSPGraph ();

    VSPGraphNode* new_node ();
    VSPGraphEdge* new_edge (VSPGraphNode* from, VSPGraphNode* to);

    void print ();
};




class VSP {
  public:
    int m_nr_elements;  // # elements
    int* m_1s;  // positive sequence
    int* m_2s;
    int* m_3s;  // negative sequence
    int* m_1s_loc;  // where is block i in ps?
    int* m_2s_loc;
    int* m_3s_loc;
    std::vector<VSPBlock*> m_blocks;  // block
    std::vector<VSPNet*> m_nets;
    double m_chip_width;
    double m_chip_depth;
    double m_chip_height;

    VSPGraph* m_XCG;
    VSPGraph* m_YCG;
    VSPGraph* m_ZCG;


    // load/store
    std::vector<int> m_1s_best;  // best solution (pos seq)
    std::vector<int> m_2s_best;  // best solution (neg seq)
    std::vector<int> m_3s_best;
    std::vector<bool> m_rotated_best;

  public:
    VSP ();
    ~VSP ();

    void clear_seq ();
    void clear_blocks ();
    void clear_CG ();
    void init (int nr_elements);
    void randomize (int seq = -1, int nr_rands = 0);
    void print (int seq = -1);
    void print_block_locations ();
    void gnu_plot (const char* file_name);
    double get_area ();
    void store_current_sequence ();
    void load_best_sequence ();
    void perturb (int type, int a, int b);
    double compute_wire();
    double compute_MIV();
    void check_fsbl(); // for checking overlapping of block coordinates

    bool create_CG ();
    void print_CG ();
    void compute_blk_locations1 ();

    // O(n^2)
   // void compute_blk_locations2 ();  // O(n log n)
};

#endif

