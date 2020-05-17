#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <map>
#include <string>
#include "st.h"

double V_blocks=0;
double W_Prv;
double D_Prv;
double H_Prv;

//std::map<char*, VSPBlock*> MapOfBlocks;

std::map<std::string, VSPBlock*> MapOfBlocks;

std::map<int, VSPBlock*> MapOfEdge;

VSPBlock::VSPBlock (int id) {

    m_id = id;
    m_name = " ";
    m_width = 0.0;
    m_depth = 0.0;
    m_height = 0.0;
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0,0;
    m_rotated = false;
};

VSPBlock::~VSPBlock () {
}

VSPNet::VSPNet(){
}
VSPNet::~VSPNet(){
}

VSPGraphNode::VSPGraphNode (int id) {
    m_id = id;

    m_edges_parents.clear ();
    m_edges_children.clear ();

    m_status = VSP_NODE_UNVISITED;
    m_data1 = NULL;
    m_loc = 0.0;
}

VSPGraphNode::~VSPGraphNode () {
}



VSPGraphEdge::VSPGraphEdge (VSPGraphNode* from, VSPGraphNode* to) {
    m_from = from;
    m_to = to;
}

VSPGraphEdge::~VSPGraphEdge () {
}



VSPGraph::VSPGraph () {
    m_nodes.clear ();
    m_edges.clear ();
    m_src_node = NULL;
    m_tail_node = NULL;
}

VSPGraph::~VSPGraph () {
    for ( std::vector<VSPGraphNode*>::iterator it = m_nodes.begin() ; it != m_nodes.end() ; ++it )
        delete (*it);

    m_nodes.clear ();

    for ( std::vector<VSPGraphEdge*>::iterator it = m_edges.begin() ; it != m_edges.end() ; ++it )
        delete (*it);

    m_edges.clear ();
}

VSPGraphNode* VSPGraph::new_node () {
    int id = int(m_nodes.size ());

    VSPGraphNode* n = new VSPGraphNode (id);
    m_nodes.push_back (n);

    return n;
}

VSPGraphEdge* VSPGraph::new_edge (VSPGraphNode* from, VSPGraphNode* to) {
    VSPGraphEdge* e = new VSPGraphEdge (from, to);
    m_edges.push_back (e);
    from->m_edges_children.push_back (e);
    to->m_edges_parents.push_back (e);

    return e;
}

void VSPGraph::print () {
    int nr_nodes = int(m_nodes.size());

    for ( int i = 0 ; i < nr_nodes ; i++ ) {
        VSPGraphNode* n;

        if ( i == (nr_nodes-2) ) {
          n = m_src_node;
          printf ("Node: src\n");
        }
        else if ( i < (nr_nodes-2) ) {
          n = m_nodes[i];
          printf ("Node: %d\n", i);
        }
        else {
          n = m_tail_node;
          printf ("Node: tail\n");
        }

        int nr_from = int(n->m_edges_parents.size());

        printf ("  From\n");

        for ( int j = 0 ; j < nr_from ; j++ ) {
          VSPGraphEdge* e = n->m_edges_parents[j];

          if ( e->m_from->m_id < int(m_nodes.size()-2) )
            printf ("    %d -> ", e->m_from->m_id);
          else if ( e->m_from->m_id == int(m_nodes.size()-2) )
            printf ("    src -> ");
          else
            printf ("    tail -> ");

          if ( e->m_to->m_id < int(m_nodes.size()-2) )
            printf ("%d\n", e->m_to->m_id);
          else if ( e->m_to->m_id == int(m_nodes.size()-2) )
            printf ("src\n");
          else
            printf ("tail\n");
        }

        int nr_to = int(n->m_edges_children.size());

        printf ("  To\n");

        for ( int j = 0 ; j < nr_to ; j++ ) {
          VSPGraphEdge* e = n->m_edges_children[j];

          if ( e->m_from->m_id < int(m_nodes.size()-2) )
            printf ("    %d -> ", e->m_from->m_id);
          else if ( e->m_from->m_id == int(m_nodes.size()-2) )
            printf ("    src -> ");
          else
            printf ("    tail -> ");

          if ( e->m_to->m_id < int(m_nodes.size()-2) )
            printf ("%d\n", e->m_to->m_id);
          else if ( e->m_to->m_id == int(m_nodes.size()-2) )
            printf ("src\n");
          else
            printf ("tail\n");
        }
    }
}



VSP::VSP () {
    m_1s = NULL;
    m_2s = NULL;
    m_3s = NULL;
    m_1s_loc = NULL;
    m_2s_loc = NULL;
    m_3s_loc = NULL;
    m_XCG = NULL;
    m_YCG = NULL;
    m_ZCG = NULL;
    m_blocks.clear ();
    m_chip_width = 0.0;
    m_chip_depth = 0.0;
    m_chip_height = 0.0;
}

VSP::~VSP () {
    clear_seq ();
    clear_blocks ();
    clear_CG ();
}

void VSP::clear_seq () {
    if ( m_1s ) {
        delete[] m_1s;
        m_1s = NULL;
      }

    if ( m_2s ) {
        delete[] m_2s;
        m_2s = NULL;
      }

    if ( m_3s ) {
        delete[] m_3s;
        m_3s = NULL;
      }

    if ( m_1s_loc ) {
        delete[] m_1s_loc;
        m_1s_loc = NULL;
      }

    if ( m_2s_loc ) {
        delete[] m_2s_loc;
        m_2s_loc = NULL;
      }

    if ( m_3s_loc ) {
        delete[] m_3s_loc;
        m_3s_loc = NULL;
      }
}

void VSP::clear_blocks () {
    for ( std::vector<VSPBlock*>::iterator it = m_blocks.begin() ; it != m_blocks.end() ; ++it )
        delete (*it);

    m_blocks.clear ();
}

void VSP::clear_CG () {
    if ( m_XCG ) {
        delete m_XCG;
        m_XCG = NULL;
      }

    if ( m_YCG ) {
        delete m_YCG;
        m_YCG = NULL;
      }

    if ( m_ZCG ) {
        delete m_ZCG;
        m_ZCG = NULL;
      }
}

void VSP::init (int nr_elements) {
    clear_seq ();
    clear_blocks ();
    clear_CG ();

    if ( nr_elements <= 0 )
        return;

    m_1s = new int[nr_elements];
    m_2s = new int[nr_elements];
    m_3s = new int[nr_elements];
    m_1s_loc = new int[nr_elements];
    m_2s_loc = new int[nr_elements];
    m_3s_loc = new int[nr_elements];
    m_nr_elements = nr_elements;

  // initialize (0 1 2 ...) (0 1 2 ...)(0 1 2 ...)
    for ( int i = 0 ; i < m_nr_elements ; i++ ) {
        m_1s[i] = i;
        m_2s[i] = i;
        m_1s_loc[i] = i;
        m_2s_loc[i] = i;

        if(i< m_nr_elements/2){
            m_3s[i] = i;
            m_3s_loc[i] = i;}
        else {
            m_3s[i] = 1.5*m_nr_elements - i -1;
            m_3s_loc[m_3s[i]] = i;
            }

      }
}

void VSP::randomize (int seq, int nr_rands) {
    int* sp;
    int* sp_loc;

      if ( seq == 0 ) {
        sp = m_1s;
        sp_loc = m_1s_loc;
      }
      else if ( seq == 1 ) {
        sp = m_2s;
        sp_loc = m_2s_loc;
      }
     else if ( seq == 2 ) {
        sp = m_3s;
        sp_loc = m_3s_loc;
      }

      if ( nr_rands <= 0 )
        nr_rands = 100*m_nr_elements;

      int a, b, t;

      if ( (seq == 0) || (seq == 1)|| (seq == 2) ) {
        for ( int i = 0 ; i < nr_rands ; i++ ) {
          a = rand() % m_nr_elements;
          b = rand() % m_nr_elements;

          t = sp[a];
          sp[a] = sp[b];
          sp[b] = t;

          sp_loc[sp[a]] = a;
          sp_loc[sp[b]] = b;
        }
      }
      else {
        for ( int i = 0 ; i < nr_rands ; i++ ) {
          a = rand() % m_nr_elements;
          b = rand() % m_nr_elements;

          t = m_1s[a];
          m_1s[a] = m_1s[b];
          m_1s[b] = t;

          m_1s_loc[m_1s[a]] = a;
          m_1s_loc[m_1s[b]] = b;

          a = rand() % m_nr_elements;
          b = rand() % m_nr_elements;

          t = m_2s[a];
          m_2s[a] = m_2s[b];
          m_2s[b] = t;

          m_2s_loc[m_2s[a]] = a;
          m_2s_loc[m_2s[b]] = b;

          a = rand() % m_nr_elements;
          b = rand() % m_nr_elements;

          t = m_3s[a];
          m_3s[a] = m_3s[b];
          m_3s[b] = t;

          m_3s_loc[m_3s[a]] = a;
          m_3s_loc[m_3s[b]] = b;
        }
      }
}


void VSP::print (int seq) {
    if ( seq != VSP_NEG ) {
        printf ("POS\n");

    for ( int i = 0 ; i < m_nr_elements ; i++ )
        printf (" %d", m_1s[i]);

        printf ("\n");
      }

    if ( seq != VSP_POS ) {
        printf ("NEG\n");

        for ( int i = 0 ; i < m_nr_elements ; i++ )
            printf (" %d", m_2s[i]);

            printf ("\n");
  }

    if ( seq != VSP_NEG ) {
        printf ("POS Loc\n");

        for ( int i = 0 ; i < m_nr_elements ; i++ )
            printf (" %d", m_1s_loc[i]);

            printf ("\n");
      }

    if ( seq != VSP_NEG ) {
        printf ("NEG Loc\n");

        for ( int i = 0 ; i < m_nr_elements ; i++ )
            printf (" %d", m_2s_loc[i]);

            printf ("\n");
  }
}


void VSP::print_block_locations () {
    for ( std::vector<VSPBlock*>::iterator it = m_blocks.begin() ; it != m_blocks.end() ; ++it ) {
        VSPBlock* b = (*it);
        printf ("  (%g , %g, %g)\n", b->m_x, b->m_y, b->m_z);
      }

      printf ("  Outline (%g , %g, %g)\n", m_chip_width, m_chip_depth, m_chip_height);
    }


void VSP::gnu_plot (const char* file_name) {
    if ( !file_name )
        return;

    FILE* fp = fopen (file_name, "w");

    if ( !fp )
        return;


   //  int nr_modules=m_modules.size();

     fprintf(fp, "plotcube([%g %g %g],[ 0 0 0], .1 ,[1 1 .7]);\n",  m_chip_width,
                   m_chip_depth,  m_chip_height);


      for (int i= 0 ; i< m_nr_elements; i++){

            float c[4];

            for ( int j = 0 ; j < 4 ; j++ )
            c[j] = ((double) rand ()) /  RAND_MAX ;



            fprintf(fp, "plotcube([%g %g %g],[ %g %g %g], %f ,[%f %f %f]);\n",   m_blocks[i]->m_width,
                   m_blocks[i]->m_depth,  m_blocks[i]->m_height,
                   m_blocks[i]->m_x, m_blocks[i]->m_y, m_blocks[i]->m_z, c[0], c[1], c[2], c[3]
                    );


  }



    fclose (fp);


}

double VSP::get_area () {
  return (m_chip_width * m_chip_depth * m_chip_height);
}

void VSP::store_current_sequence () {
    m_1s_best.clear ();
    m_2s_best.clear ();
    m_3s_best.clear ();

    m_rotated_best.clear ();

    for ( int i = 0 ; i < m_nr_elements ; i++ ) {
        m_1s_best.push_back (m_1s[i]);
        m_2s_best.push_back (m_2s[i]);
        m_3s_best.push_back (m_3s[i]);
        m_rotated_best.push_back (m_blocks[i]->m_rotated);
  }
}

void VSP::load_best_sequence () {
    if ( (int(m_1s_best.size()) != m_nr_elements) || (int(m_2s_best.size()) != m_nr_elements) ||(int(m_3s_best.size()) != m_nr_elements) || (int(m_rotated_best.size()) != m_nr_elements) )
        return;

    for ( int i = 0 ; i < m_nr_elements ; i++ ) {
        m_1s[i] = m_1s_best[i];
        m_2s[i] = m_2s_best[i];
        m_3s[i] = m_3s_best[i];

        m_blocks[i]->m_rotated = m_rotated_best[i];
      }

 // printf("Loading Best Seq\n");
}

void VSP::perturb (int type, int a, int b) {
    if ( type == 0 ) {  // type-1: swap two cells in the positive sequence
        int t = m_1s[a];
        m_1s[a] = m_1s[b];
        m_1s[b] = t;

        m_1s_loc[m_1s[a]] = a;
        m_1s_loc[m_1s[b]] = b;
      }
    else if ( type == 1 ) {  // type-2: swap two cells in the negative sequence
        int t = m_2s[a];
        m_2s[a] = m_2s[b];
        m_2s[b] = t;

        m_2s_loc[m_2s[a]] = a;
        m_2s_loc[m_2s[b]] = b;
      }
    else if ( type == 2 ) {  // type-3: swap two cells in the  sequence 3
        int t = m_3s[a];
        m_3s[a] = m_3s[b];
        m_3s[b] = t;

        m_3s_loc[m_3s[a]] = a;
        m_3s_loc[m_3s[b]] = b;
      }

    else if ( type == 3 ) {  // type-4: swap two cells in the all sequences
        int t = m_1s[a];
        m_1s[a] = m_1s[b];
        m_1s[b] = t;

        m_1s_loc[m_1s[a]] = a;
        m_1s_loc[m_1s[b]] = b;

        t = m_2s[a];
        m_2s[a] = m_2s[b];
        m_2s[b] = t;

        m_2s_loc[m_2s[a]] = a;
        m_2s_loc[m_2s[b]] = b;

        t = m_3s[a];
        m_3s[a] = m_3s[b];
        m_3s[b] = t;

        m_3s_loc[m_3s[a]] = a;
        m_3s_loc[m_3s[b]] = b;

      }
      else if ( type == 4) {  // type-5: rotate a block (block a)

        int rnd_h= rand() % n_Layers + 1;

        double  AR_1 = rand()% 700+600 ;
        double AR = AR_1/1000;  // AR ratio from .7 to 1.3


        W_Prv = m_blocks[a]->m_width;
        D_Prv = m_blocks[a]->m_depth ;
        H_Prv = m_blocks[a]->m_height;

        double V_m = W_Prv * D_Prv * H_Prv;
        m_blocks[a]->m_height= rnd_h ;
        m_blocks[a]->m_width = sqrt((V_m * AR /rnd_h)) ;
        m_blocks[a]->m_depth = m_blocks[a]->m_width /AR ;


   // m_blocks[a]->m_rotated = !(m_blocks[a]->m_rotated);
  }
}

bool VSP::create_CG () {
    clear_CG ();

    m_XCG = new VSPGraph ();
    m_YCG = new VSPGraph ();
    m_ZCG = new VSPGraph ();


  // internal nodes
    for ( int i = 0 ; i < m_nr_elements ; i++ ) {
        m_XCG->new_node ();
        m_YCG->new_node ();
        m_ZCG->new_node ();

      }

  // source nodes
    m_XCG->m_src_node = m_XCG->new_node ();
    m_YCG->m_src_node = m_YCG->new_node ();
    m_ZCG->m_src_node = m_ZCG->new_node ();


  // tail nodes
    m_XCG->m_tail_node = m_XCG->new_node ();
    m_YCG->m_tail_node = m_YCG->new_node ();
    m_ZCG->m_tail_node = m_ZCG->new_node ();

  // node index (block id -> index in the sequence pair)
    std::vector<int> node_index_in_1s (m_nr_elements);
    std::vector<int> node_index_in_2s (m_nr_elements);
    std::vector<int> node_index_in_3s (m_nr_elements);

    for ( int i = 0 ; i < m_nr_elements ; i++ ) {
        node_index_in_1s[m_1s[i]] = i;
        node_index_in_2s[m_2s[i]] = i;
        node_index_in_3s[m_3s[i]] = i;

        m_XCG->m_nodes[i]->m_data1 = (void*) m_blocks[i];
        m_YCG->m_nodes[i]->m_data1 = (void*) m_blocks[i];
        m_ZCG->m_nodes[i]->m_data1 = (void*) m_blocks[i];

  }

  // pair between block i and block k
    int x1, x2, x3, y1, y2, y3;

    for ( int i = 0 ; i < m_nr_elements ; i++ ) {
        for ( int k = i+1 ; k < m_nr_elements ; k++ ) {
          x1 = node_index_in_1s[i];
          x2 = node_index_in_2s[i];
          x3 = node_index_in_3s[i];

          y1 = node_index_in_1s[k];
          y2 = node_index_in_2s[k];
          y3 = node_index_in_3s[k];


          if ( (x2 < y2) && (x3 > y3) )       //(..i..k../..k..i..) (..i..k..) (..k..i..) => k is left to i
            m_XCG->new_edge (m_XCG->m_nodes[k], m_XCG->m_nodes[i]);
          else if ( (x2 > y2) && (x3 < y3) )  // (..i..k../..k..i..)(..k..i..) (..i..k..) => k is right to i
            m_XCG->new_edge (m_XCG->m_nodes[i], m_XCG->m_nodes[k]);
          else if ( (x1 > y1) && (x2 > y2) && (x3 > y3) )  // (..i..k..)(..i..k..) (..i..k..)  => k is front i
            m_YCG->new_edge (m_YCG->m_nodes[k], m_YCG->m_nodes[i]);
          else if ( (x1 < y1) && (x2 < y2) && (x3 < y3) )  // (..k..i..)(..k..i..) (..k..i..) => K is rear of i
            m_YCG->new_edge (m_YCG->m_nodes[i], m_YCG->m_nodes[k]);
          else if ( (x1 > y1) && (x2 < y2) && (x3 < y3) )  // (..k..i..)(..i..k..) (..i..k..)  => k is above of i
            m_ZCG->new_edge (m_ZCG->m_nodes[i], m_ZCG->m_nodes[k]);
          else if ( (x1 < y1) && (x2 > y2) && (x3 > y3) )  // (..i..k..)(..k..i..) (..k..i..) => K is below of i
            m_ZCG->new_edge (m_ZCG->m_nodes[k], m_ZCG->m_nodes[i]);
          else
            printf ("Error!!!!!!\n");
        }
  }

  // source and tail
      for ( int i = 0 ; i < m_nr_elements ; i++ ) {
        if ( m_XCG->m_nodes[i]->m_edges_parents.size() == 0 )
          m_XCG->new_edge (m_XCG->m_src_node, m_XCG->m_nodes[i]);
        if ( m_XCG->m_nodes[i]->m_edges_children.size() == 0 )
          m_XCG->new_edge (m_XCG->m_nodes[i], m_XCG->m_tail_node);

        if ( m_YCG->m_nodes[i]->m_edges_parents.size() == 0 )
          m_YCG->new_edge (m_YCG->m_src_node, m_YCG->m_nodes[i]);
        if ( m_YCG->m_nodes[i]->m_edges_children.size() == 0 )
          m_YCG->new_edge (m_YCG->m_nodes[i], m_YCG->m_tail_node);

        if ( m_ZCG->m_nodes[i]->m_edges_parents.size() == 0 )
          m_ZCG->new_edge (m_ZCG->m_src_node, m_ZCG->m_nodes[i]);
        if ( m_ZCG->m_nodes[i]->m_edges_children.size() == 0 )
          m_ZCG->new_edge (m_ZCG->m_nodes[i], m_ZCG->m_tail_node);
      }

  return true;
}

void VSP::print_CG () {
    for ( int i = 0 ; i < 3 ; i++ ) {
        VSPGraph* g;

        if ( i == 0 ) {
            g = m_XCG;
          printf ("XCG\n");
        }
        else if ( i == 1) {
          g = m_YCG;
          printf ("YCG\n");
        }
        else {
          g = m_ZCG;
          printf ("ZCG\n");
        }

        for ( std::vector<VSPGraphNode*>::iterator it = g->m_nodes.begin() ; it != g->m_nodes.end() ; ++it ) {
            VSPGraphNode* n = (*it);

            printf ("  Node %d\n", n->m_id);

        for ( std::vector<VSPGraphEdge*>::iterator ite = n->m_edges_parents.begin() ; ite != n->m_edges_parents.end() ; ++ite ) {
            VSPGraphEdge* e = (*ite);

            printf ("    e %d -> %d\n", e->m_from->m_id, e->m_to->m_id);
          }

        for ( std::vector<VSPGraphEdge*>::iterator ite = n->m_edges_children.begin() ; ite != n->m_edges_children.end() ; ++ite ) {
            VSPGraphEdge* e = (*ite);

            printf ("    e %d -> %d\n", e->m_from->m_id, e->m_to->m_id);
      }
    }
  }
}

static void sp_visit_parents (VSPGraph* graph, VSPGraphNode* cur_node, int XYZ) {
    double max_loc = 0.0;

    for ( std::vector<VSPGraphEdge*>::iterator ite = cur_node->m_edges_parents.begin() ; ite != cur_node->m_edges_parents.end() ; ++ite ) {
        VSPGraphEdge* edge = (*ite);
        VSPGraphNode* parent = edge->m_from;

        if ( parent == graph->m_src_node )
        continue;

        if ( parent->m_status == VSP_NODE_UNVISITED )
          sp_visit_parents (graph, parent, XYZ);

        double loc;

        VSPBlock* p = (VSPBlock*) (parent->m_data1);

        if ( XYZ==0 ) {

            loc = parent->m_loc + ((VSPBlock*) parent->m_data1)->m_width;
        }
        else if ( XYZ == 1){

            loc = parent->m_loc + ((VSPBlock*) parent->m_data1)->m_depth;
        }
        else {
            loc = parent->m_loc + ((VSPBlock*) parent->m_data1)->m_height;
        }

        if ( loc > max_loc )
          max_loc = loc;
      }

      cur_node->m_status = VSP_NODE_VISITED;

      if ( max_loc > cur_node->m_loc )
        cur_node->m_loc = max_loc;

      if ( cur_node != graph->m_tail_node ) {
        if ( XYZ==0 )
          ((VSPBlock*) cur_node->m_data1)->m_x = max_loc;
        else if ( XYZ==1 )
          ((VSPBlock*) cur_node->m_data1)->m_y = max_loc;
        else if ( XYZ == 2)
          ((VSPBlock*) cur_node->m_data1)->m_z = max_loc;
      }
}

void VSP::compute_blk_locations1 () {
      create_CG ();
//  print ();
//  print_CG ();

      sp_visit_parents (m_XCG, m_XCG->m_tail_node, 0);
      sp_visit_parents (m_YCG, m_YCG->m_tail_node, 1);
      sp_visit_parents (m_ZCG, m_ZCG->m_tail_node, 2);

      m_chip_width = m_XCG->m_tail_node->m_loc;
      m_chip_depth = m_YCG->m_tail_node->m_loc;
      m_chip_height = m_ZCG->m_tail_node->m_loc;
}


void VSP::check_fsbl(){


    for(int i= 0; i< m_nr_elements ;i++){

        float x3 = m_blocks[i]->m_x;
        float y3 = m_blocks[i]->m_y;
        float z3 = m_blocks[i]->m_z;

        float x4 = m_blocks[i]->m_x + m_blocks[i]->m_width;
        float y4 = m_blocks[i]->m_y + m_blocks[i]->m_depth;
        float z4 = m_blocks[i]->m_z + m_blocks[i]->m_height;


        for(int j=i+1; j< m_nr_elements; j++){

            float x1 = 0;
            float x2 = 0;
            float y1 = 0;
            float y2 = 0;
            float z1 = 0;
            float z2 = 0;

            x1 = m_blocks[j]->m_x;
            x2 = m_blocks[j]->m_x + m_blocks[j]->m_width;

            y1 = m_blocks[j]->m_y;
            y2 = m_blocks[j]->m_y + m_blocks[j]->m_depth;

            z1 = m_blocks[j]->m_z;
            z2 = m_blocks[j]->m_z + m_blocks[j]->m_height;

            int overlap_x = 0;
            int overlap_y = 0;
            int overlap_z = 0;

            if (((x1<x3)&(x2>x3))|((x1<x4)&(x2>x4))|((x1>x3)&(x1<x4))|((x2>x3)&(x2<x4))){

                 overlap_x = 1;
            }

            if (((y1<y3)&(y2>y3))|((y1<y4)&(y2>y4))|((y1>y3)&(y1<y4))|((y2>y3)&(y2<y4))){

                 overlap_y = 1;
            }

            if (((z1<z3)&(z2>z3))|((z1<z4)&(z2>z4))|((z1>z3)&(z1<z4))|((z2>z3)&(z2<z4))){

                 overlap_z = 1;
            }

            if ((x1==x4)|(x2==x3)){

                overlap_x =0;
            }

            if ((y1==y4)|(y2==y3)){

                overlap_y =0;
            }
            if ((z1==z4)|(z2==z3)){

                overlap_z =0;
            }


            if ((overlap_x + overlap_y + overlap_z)==3){
                printf("\n Error : Overlap i: %d j:%d\n", i,j);
                printf("\n (%g, %g, %g) (%g, %g, %g) ",x3,y3,z3,x4,y4,z4);
                printf("\n (%g, %g, %g) (%g, %g, %g) \n", x1,y1,z1,x2,y2,z2);

            }
        }
    }
}

double VSP::compute_wire(){


    double total_length=0;


    for(int i=0;i< m_nets.size();i++){


        if( m_nets[i]->m_net.size()>1){



            double min_x = m_nets[i]->m_net[0]->m_x;
            double min_y = m_nets[i]->m_net[0]->m_y;

            double max_x = m_nets[i]->m_net[0]->m_x + m_nets[i]->m_net[0]->m_width;
            double max_y = m_nets[i]->m_net[0]->m_y + m_nets[i]->m_net[0]->m_depth;

            double x1 = 0;
            double x2 = 0;
            double y1 = 0;
            double y2 = 0;


            for(int j=1; j< m_nets[i]->m_net.size(); j++)
            {

                x1 = m_nets[i]->m_net[j]->m_x;
                x2 = m_nets[i]->m_net[0]->m_x + m_nets[i]->m_net[j]->m_width;
                y1 = m_nets[i]->m_net[j]->m_y;
                y2 = m_nets[i]->m_net[0]->m_y + m_nets[i]->m_net[j]->m_depth;

                if(x1<min_x){
                    min_x = x1;
                }

                if(x2>max_x){
                    max_x = x2;
                }

                if(y1<min_y){
                    min_y = y1;
                }

                if(y2>max_y){
                    max_y = y2;
                }

            }

            total_length = total_length + (max_x - min_x + max_y - min_y)/2;

        }

    }

    return total_length;
}


double VSP::compute_MIV(){

    double total_MIV=0;
    int n1=0;
    int n2=0;
    int n3=0;
    int n4=0;

    for(int i=0;i< m_nets.size();i++){


        if( m_nets[i]->m_net.size()>1){

            double max_z = m_nets[i]->m_net[0]->m_z;

            double min_z = m_nets[i]->m_net[0]->m_z + m_nets[i]->m_net[0]->m_height;

            double z1 = 0;
            double z2 = 0;

            for(int j=1; j< m_nets[i]->m_net.size(); j++){

                z1 = m_nets[i]->m_net[j]->m_z;
                z2 = m_nets[i]->m_net[j]->m_z + m_nets[i]->m_net[j]->m_height;


                if(z2<min_z){
                    min_z = z2;
                }

                if(z1>max_z){
                    max_z = z1;
                }

            }

            double miv = (max_z - min_z);

            if (miv<0){
                miv=0;
           }

            total_MIV = total_MIV + miv ;

        }

    }

    return total_MIV;

}



bool open_blk (VSP& sp, const char* file_name) {

      FILE* fp = fopen (file_name, "r");

      if ( !fp )
        return false;

      char szBuf[1024];
      memset (szBuf, 0x00, sizeof(szBuf));

      if ( fgets (szBuf, sizeof(szBuf)-1, fp) == NULL ) {
        fclose (fp);
        return -1;
      }


      int nr_blks = atoi (szBuf);
      const char* tok = " \t\r\n";
      bool err = false;

      sp.init (nr_blks);




      for ( int i = 0 ; i < nr_blks ; i++ ) {

        memset (szBuf, 0x00, sizeof(szBuf));
        if ( fgets (szBuf, sizeof(szBuf)-1, fp) == NULL ) {

            break;
        }

        char* pch = strtok (szBuf, tok);

        if(!pch) // if nothing is in 'pch(block)' (because there is some empty line in given file)
			continue;



		std::string bname = pch;// it chould be a block name
		//printf("%s \n", pch);

		pch = strtok(NULL, tok);

		if(!pch)
			continue;

		if(strcmp(pch, "softrectangular") != 0)
			continue; // go to beginning of the 'while' loop, not go down.


		pch= strtok(NULL, tok );
		strtok(NULL, tok );
		strtok(NULL, tok );

        int rnd_h= rand() % n_Layers + 1;
        double h = rnd_h;
        double w = sqrt((atoi(pch) / h));
        double d = w;

        double V = w*d*h;
        V_blocks = V_blocks + V;

        if ( (w <= 0.0) || (d <= 0.0)||(h <= 0.0) ) {

          break;

        }

        VSPBlock* temp_module;
        temp_module = new VSPBlock(i);

        temp_module->m_width=w;
        temp_module->m_height=h;
        temp_module->m_depth=d;
        temp_module->m_name=bname;
        //p,w,d,h,bname);

        (sp.m_blocks).push_back(temp_module);


        MapOfBlocks.insert(std::make_pair(bname, temp_module));

}

  if ( err ) {
    fclose (fp);
    return -1;
  }


  fclose (fp);

  return true;
}


bool open_net (VSP& sp, const char* file_name) {

    FILE* np;
    np = fopen (file_name, "r");
    if ( np == NULL ) {
    printf ("cannot open the file\n");
    return -3;
    }


    int x=1;

while(1)
	{   char buf[16];
		memset(buf, 0x00, sizeof(buf));

		if(fgets(buf, sizeof(buf)-1, np) == NULL)
		{
			break;

		}

		const char* tok = " \t\r\n";
        bool err = false;
		char* pch = strtok(buf, tok);

		if(!pch)
			continue;
        //puts(pch);


		if(strcmp(pch, "NetDegree") != 0)
			continue;

            //puts(pch);
        pch=strtok(NULL, ":" );
       // printf("%s\n",pch );

		int n_net= atoi(pch);
		//printf("%d \n",n_net);

        VSPNet* temp_net=new VSPNet();


		for(int i=0;i<n_net;i++)
        {

        //printf("%d \n", i);
		memset(buf, 0x00, sizeof(buf));

		if(fgets(buf, sizeof(buf)-1, np) == NULL)
		{
			break;

		}
        //puts(buf);
		pch=strtok(buf,tok);

		VSPBlock* temp_block;
        temp_block=MapOfBlocks.find(pch)->second ;
        //std::cout<< temp<<"\n";


        if(temp_block != NULL)
        {
            (temp_net->m_net).push_back(temp_block);
            temp_net->m_num=x;
            (temp_block->conn_net).push_back(temp_net);
            //std::string str = std::to_string(x);
            //temp_cnet->m_name = "net"+ std::to_string((long double)x);           //std::cout<<temp_cnet->m_net<<"\n";
        }
        else continue;

		//std:: map<std::string, CBlock*>::iterator itr;
		//itr=MapOfBlocks.find(pch);

        }

      //  std::cout<<"net# "<<x<<"    "<<temp_cnet<<"\n";

        x++;

        (sp.m_nets).push_back(temp_net);
	}

  fclose (np);

  return true;
}


static void get_two_random_numbers (int& a, int& b, int max) {
  a = rand() % max;

  while ( 1 ) {
    b = rand() % max;

    if ( a != b )
      break;
  }
}



static void metropolis (VSP& sp, double T, int M, double& best_cost, double AR_min, double AR_max) {
  // 1. Declare five floating-point variables (double cur_cost, new_cost, delta_h, r, AR)
    float cur_cost; // cur_cost: the cost of the current solution
    float new_cost;//
    float delta_h;// delta_h: cost(NewS) - cost(S)
    float r = rand ()% 1;// r: a random number between (0, 1)
    float AR;// AR: aspect ratio of the new (perturbed) solution
    float alpha = 1;
    float beta = 20;

  // 2*. Perturb the current solution #M times.
    for ( int i = 0 ; i < M ; i++ ) {
    // 3*. Get the current cost (chip area)
        sp.compute_blk_locations1 ();

        cur_cost = alpha*sp.get_area () + beta*sp.compute_wire();

        int move = rand() % 5;  // five different types (0, 1, 2, 3, 4)
        int a, b;
        get_two_random_numbers (a, b, sp.m_nr_elements);  // get two random numbers in [0, # blocks - 1]

        sp.perturb (move, a, b);  // perturb the current solution.


        sp.compute_blk_locations1 ();

        new_cost = alpha*sp.get_area () + beta*sp.compute_wire();

        delta_h = new_cost-cur_cost; // 6. Compute delta_h (= cost(NewS) - cost(S)).

        AR=sp.m_chip_width/sp.m_chip_depth;

        // 8*. Generate a random number in (0, 1) and save it into r.
        r = ((double) rand ()) / RAND_MAX ;

        if ( ((delta_h < 0) || (r < exp(-1.0*delta_h/T))) && ((AR_min <= AR) && (AR <= AR_max) && (sp.m_chip_height <= n_Layers)) ) {  // accept
          // 10*. Check whether this is the best solution.

          if ( new_cost < best_cost ) {
                best_cost = new_cost;
                    printf("best area: %g   best wire_length: %f   MIV: %g \n", sp.get_area (), sp.compute_wire() , sp.compute_MIV());
                    sp.check_fsbl();
                    sp.store_current_sequence();
                    //sp.gnu_plot ("final.plt");
              }
        }

        else {  // reject
            if ( move == 4){
                sp.m_blocks[a]->m_height = H_Prv;
                sp.m_blocks[a]->m_width = W_Prv;
                sp.m_blocks[a]->m_depth = D_Prv;
            }
            else {
                sp.perturb(move,a,b);
            }



    }
  }
}






static void simulated_annealing (VSP& sp) {


    sp.store_current_sequence ();

    double T0=100;

    double zeta=.99;

    double theta=.99;

    double AR_min=.7; double AR_max=1.3;

    double T=T0;

    int M=500*sp.m_nr_elements* log(1+ n_Layers);

    int my_time=0;

    int max_time=500000*sp.m_nr_elements* log(1+ n_Layers);

    double best_cost=100000000;


    while ( (my_time <= max_time) && (T >= 1.0) ) {
        printf ("T: %g\n", T);  // always print out the current temperature to see the progress.
        metropolis (sp, T, M, best_cost, AR_min, AR_max);  // call metropolis.
        my_time=my_time+M; // 11. Increase my_time by M. (Time = Time + M;)

        T=zeta*T;// 12. Reduce the temperature by alpha. (T = alpha * T; where alpha = cooling rate (alpha < 1))

        M=theta*M; // 13. Reduce the # moves to try. (M = beta * M; where beta < 1)
      }


  sp.load_best_sequence ();


  sp.compute_blk_locations1 ();
}





int main (int argc, char* argv[]) {
      if ( argc != 3 ) {
        printf ("Usage: main  <input>\n");
        return -1;
      }

      VSP sp;

      open_blk (sp, argv[1]);

      printf("\n sum of volumes: %g\n", V_blocks);
      open_net (sp, argv[2]);

      // this initializes the sequence pair
      srand (time(NULL));

      sp.compute_blk_locations1 ();
      double init_area = sp.get_area ();
      printf ("init cost:  %g\n", init_area);
      printf ("init height: %g\n", sp.m_chip_height);  // get the initial floorplan area.
      sp.gnu_plot ("init.plt");  // dump the current floorplan into a file.

      simulated_annealing (sp);
      sp.gnu_plot ("final.plt");


      return 0;
}

