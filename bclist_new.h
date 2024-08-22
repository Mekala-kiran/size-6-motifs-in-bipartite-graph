#ifndef bclist_new_h
#define bclist_new_h

#include <iostream>
#include<fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <time.h>
#include <cmath>
#include <random>
#include <sys/resource.h>
#include <tuple> 

#define COUNT_ONLY

#define USE_CORE_REDUCTION

#define SZ(x) ((int)x.size())
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

#define NTWOHOPS 10000

#define SAMPLE_RATE 0.01

#define MAX_N 10000000
#define MAX_K 20
#define IS_DEBUGGING

using namespace std;

typedef int long long lint;
#define ll long long
double time_counter0 = 0.0;
double time_counter1 = 0.0;

lint combination_cache[MAX_N+1][MAX_K+1] = {0};
lint my_factorial_cache[MAX_N+1][MAX_K+1] = {0};

class VertexDegree{
public:
    int vertex;
    int degree;
    inline bool operator < (const VertexDegree &other) const {
        return degree > other.degree || (degree == other.degree && vertex < other.vertex);
    }
    
    VertexDegree();
    VertexDegree(int v, int d);
    ~VertexDegree();
};
bool cmpByDegAsc(VertexDegree &v1, VertexDegree &v2){
    return v1.degree < v2.degree || (v1.degree == v2.degree && v1.vertex < v2.vertex);
}

class LargeBiclique{
public:
    vector<int> left_vertices;
    vector<int> right_vertices;
    int size;
    
    inline bool operator < (const LargeBiclique &other) const {
        return size > other.size;
    }
    
    LargeBiclique();
    LargeBiclique(vector<int> &l_vertices, vector<int> &r_vertices);
    ~LargeBiclique();
};

class Tools{
public:
    static vector<int> intersection(vector<int> &vec1, vector<int> &vec2);
    static vector<int> intersection(vector<int> &vec1, vector<int> &vec2, int offset1, int offset2);
    static vector<int> merge(vector<int> &vec1, vector<int> &vec2);
    static int intersection_count(vector<int> &vec1, vector<int> &vec2);
    static lint choose(lint n, lint k);
    static lint my_factorial(lint n, lint k);
    static string toString(vector<int> &vec);
    static vector<int> toVector(const string &str);
};

class SpecialBigraph {
public:
    int num_edges;
    int num_vertices;
    
    int n_vertices[2];
    map<int, int> vertices[2];

    
    vector<int> vertices_in_left;
    vector<int> vertices_in_right;
    
    vector<int> all_vertices;
    
    //unordered_map<int, int> vertexranks;
    unordered_map<int, int> vertexids;
    
    int *deg;
    vector< vector <int> > adj_vec;
    //vector< vector <int> > two_hop_adj_vec;
    int **two_hop_adj_vec;
    int *two_hop_adj_size;
    int *two_hop_adj_maxsize;

	set <vector<int>> edges;
map<pair<int,int>,int> edge_sign; 
vector < pair <int, int>> list_of_edges;
    
    vector<int> clique_vertices_in_left;
    
    vector< pair<vector<int>, vector<int> > > bcliques;

    
    lint largest_index_in_partition[2];
    
    int p;
    int q;
    int n;
    int priority;
    int use_cost_model;
    
    int *ns;
    int *d;
    int *cd;
    int *adj;
    int **two_hop_d;
    int *two_hop_cd;
    int *two_hop_adj;
    int *lab;
    int **sub;
    int **op_vertices;
    int *op_size;
    lint result_count;
    
private:
    string data_file_path;
    
    bool anchor_left;
    
    void clear_everything();
    bool all_num(string &s);
    void add_vertex(int A, int side);
    void add_edge(int &A, int &B, int &sign); // Modified to handle signed edges

    void get_index(int &A, int side);

    void reformat_graph();
    
    void collect_two_hop_adj();
    void collect_two_hop_adj_bk();
    
    void trim_graph_by_core();
    
    double estimate_cost(int side);
    
    void sort_vertices(int strategy);
    void sort_vertices_bk(int strategy);
    
    void pqclique(int l);

public:
    // tools
    void print(vector<int> vec);
    void print_adj();
    void print_adj(int i);
    void print_two_hop_adj();
    void print_edges();
    void print_deg();
    void print_bclique(ostream &os, vector<int> &left_vertices, vector<int> &right_vertices);
    void print_map(unordered_map<int, int> ids);
    void print_results();
    void print_graph_with_signs();

    
    SpecialBigraph(string in_path, int p_value, int q_value, int priority_strategy, int use_cost_model_flag=1);
    ~SpecialBigraph();
    void read_graph();
    
    bool prepare_graph();
    
    void count_butterflies();
    
    void listing_cliques();
   
};
#endif /* bclist_new_h */
