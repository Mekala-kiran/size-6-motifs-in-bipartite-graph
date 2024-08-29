#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <time.h>
#include <cstdio>
#include <cassert>
#include <cstdio>
#include <stdio.h>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <sstream>
#include <chrono>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/random_device.hpp>
#include <boost/variant.hpp>

using namespace std;
using namespace boost::random;

#define SZ(x) ((int)x.size())
#define ll long long
#define ull unsigned long long
#define ld long double
#define eps 1e-11
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

const int ITER_VER = 2200;
const ll shift = 1000 * 1000 * 1000LL;
const double TIME_LIMIT = 20;
const int N_WEDGE_ITERATIONS = 2 * 1000 * 1000 * 10;
const int ITERATIONS_SAMPLING = 5;
const int N_SPARSIFICATION_ITERATIONS = 5;
const int TIME_LIMIT_SPARSIFICATION = 10000; // !half an hour
const int N_FAST_EDGE_BFC_ITERATIONS = 2100; // used for fast edge sampling
const int N_FAST_WEDGE_ITERATIONS = 50; // used for fast wedge sampling

char output_address [2000] ;
char *input_address;


set <vector<int>> edges;
map<pair<int,int>,int> edge_sign,new_edge_sign; 
vector < pair <int, int>> list_of_edges;
map < int, int > vertices [2];
vector <int> index_map;
vector <int> vertices_in_left;
vector <int> vertices_in_right;
vector < vector <int> > adj,newadj;
vector < vector < int > > sampled_adj_list;
vector <bool> visited;
vector <int> list_of_vertices;
vector <int> vertex_counter;
unordered_map<int,int> degree_left,degree_right;

ll count_wedge;
ll n_vertices;
ll n_edges;
ld exact_n_bf;
ld exact_n_bf_signed;
ld exact_n_bf_unsigned;
ll n_wedge_in_partition[2];
ll largest_index_in_partition[2];

vector <int> clr;
vector <int> hashmap_C;
vector <ll> sum_wedges;
vector <ll> sum_deg_neighbors;
vector <int> aux_array_two_neighboorhood;

void clear_everything() {
	largest_index_in_partition[0] = largest_index_in_partition[1] = 0;
	n_vertices = 0;
	n_edges = 0;
	edges.clear();
	vertices[0].clear(); vertices[1].clear();
	index_map.clear();
	vertices_in_left.clear();
	vertices_in_right.clear();
	adj.clear();
	sampled_adj_list.clear();
	visited.clear();
	//list_of_edges.clear();
	vertex_counter.clear();
	clr.clear();
	hashmap_C.clear();
	sum_wedges.clear();
	sum_deg_neighbors.clear();
	aux_array_two_neighboorhood.clear();
}

void resize_all() {
	clr.resize(n_vertices);
	hashmap_C.resize(n_vertices);
	aux_array_two_neighboorhood.resize(n_vertices);
	sum_wedges.resize(n_vertices);
	visited.resize(n_vertices);
	index_map.resize(n_vertices);
	sum_deg_neighbors.resize(n_vertices);
}



void add_vertex(int A, int side) {
	if (vertices[side].find(A) == vertices[side].end()) {
		if (side == 0) vertices_in_left.push_back(A);
		else vertices_in_right.push_back(A);
		vertices[side][A] = 1;
	}
}

void add_edge(int &A, int &B, int &sign) {
	add_vertex(A, 0);
	add_vertex(B, 1);

	if (edges.find({A,B,sign}) == edges.end()) {
		edges.insert({A,B,sign});

		n_edges++;
	}
}

void get_index(int &A, int side) {
	if (vertices[side].find(A) == vertices[side].end()) {
		vertices[side][A] = largest_index_in_partition[side] ++ ;
		/*print A and vertices[side][A] to get the generated id for the vertex A*/
		
		
		//cout<<"i am A"<<'\t'<<A<<'\n';
		// cout<< A << " " << side << " " << vertices[side][A]<<'\n';		
	}
	A = vertices[side][A];
}

bool all_num(string &s) {
	for (int i = 0; i < SZ(s); i++) if ((s[i] >= '0' && s [i] <= '9') == false) return false;
	return true;
}

void get_graph() {
	freopen(input_address, "r", stdin); //tries to open a file with a file stream that is associated with another opened file.
	printf("[%d:] File Opened\n", __LINE__);
	string s;
	cin.clear(); //clears the error flag on cin 
	while (getline(cin, s)) { //to read a string or a line from an input stream
 		stringstream ss; ss << s; //A stringstream associates a string object with a stream allowing you to read from the string as if it were a stream (like cin). To use stringstream, we need to include sstream header file.
 		
		vector <string> vec_str; 
		for (string z; ss >> z; vec_str.push_back(z));//push elements into a vector from the back
		if (SZ(vec_str) >= 2) {
			bool is_all_num = true;
			for (int i = 0; i < min (2, SZ(vec_str)) ; i++) is_all_num &= all_num(vec_str[i]);
			if (is_all_num) {
				int A, B, sign;
				ss.clear(); ss << vec_str[0]; ss >> A;
				ss.clear(); ss << vec_str[1]; ss >> B;
                ss.clear(); ss << vec_str[2]; ss >> sign;
				add_edge(A, B, sign);
				
			}
		}
	}
	printf("[%d:] File Reading Done!\n", __LINE__);
	vertices[0].clear();
	vertices[1].clear();
	largest_index_in_partition[0] = 0;
	largest_index_in_partition[1] = SZ(vertices_in_left);
	n_vertices = SZ(vertices_in_left) + SZ(vertices_in_right);
	adj.resize(n_vertices, vector <int> ());
	for (auto edge : edges) {
		int A = edge[0];  
		int B = edge[1];
        	int sign = edge[2];
		get_index(A, 0);
		get_index(B, 1);
		adj[A].push_back(B);
		adj[B].push_back(A);
        //list_of_edges.push_back(make_pair(A, B)); 
        edge_sign[{A, B}] = sign;
		edge_sign[{B, A}] = sign;
		
		
	}
	resize_all();

	n_wedge_in_partition[0] = 0;
	for (int i = 0; i < largest_index_in_partition[0]; i++) {
		n_wedge_in_partition[0] += (((ll)SZ(adj[i])) * (SZ(adj[i]) - 1)) >> 1;
	}
	n_wedge_in_partition[1] = 0;
	for (int i = largest_index_in_partition[0]; i < largest_index_in_partition[1]; i++) {
		n_wedge_in_partition[1] += ((ll)SZ(adj[i]) * (SZ(adj[i]) - 1)) >> 1;
	}
	for (int i = 0; i < n_vertices; i++) {
		sort(adj[i].begin(), adj[i].end());
		sum_deg_neighbors[i] = 0;
		for (auto neighbor : adj[i]) {
			sum_deg_neighbors[i] += SZ(adj[neighbor]);
		}
	}
	//cerr << " for test # edges :: " << SZ(list_of_edges) << " left :: " << SZ(vertices_in_left) << " right :: "  << SZ(vertices_in_right) << endl;
	//sort(list_of_edges.begin(), list_of_edges.end());
	edges.clear();
	fclose(stdin);
}

/*This function returns 1 if priority(u) < priority(v), otherwise it returns 0*/
//int priority(int u, int v){
//}

void read_the_graph(char *fin) {
	clear_everything();
	//cerr << " Insert the input (bipartite network) file location" << endl;
	//cerr << " >>> "; 
	input_address = fin;
	//cerr << " Insert the output file" << endl;
	//cerr << " >>> "; cin >> output_address;
	//freopen(output_address, "w", stdout);
	cerr << " ---------------------------------------------------------------------------------------------------------------------- \n";
	cerr << "| * Note that edges should be separated line by line.\n\
| In each line, the first integer number is considered as a vertex in the left partition of bipartite network, \n\
| and the second integer number is a vertex in the right partition. \n\
| In addition, multiple edges are removed from the given bipartite network.\n\
| Also, note that in this version of the source code, we did NOT remove vertices with degree zero.\n";
	cerr << " ---------------------------------------------------------------------------------------------------------------------- \n";

	cerr << " Processing the graph ... (please wait) \n";

	get_graph();   //function() from 27th line

	cout << " -------------------------------------------------------------------------- \n";
	cout << "Input graph: " << input_address << "\n";
	cout << " The graph is processed - there are " << n_vertices << " vertices and " << n_edges << " edges  \n";
	cout << " -------------------------------------------------------------------------- \n";
}





/*ll all_bal_bat_counting(vector < vector <int> > &graph) {
ll res=0;
int v;
	for(int u=0; u < vertices_in_left.size(); u++){
		
		unordered_map<int, int> count_wedge_with_signs_0;
		unordered_map<int, int> count_wedge_with_signs_1; 
		unordered_map<int, int> count_wedge_with_signs_2; 
		set<int> n_w;
		for(int j = 0; j < SZ(graph[u]); j++)
		{
			 v = graph[u][j];
			int sign_sum = edge_sign[{u,v}]; 
			//if(1||SZ(graph[v]) < SZ(graph[u]) || ((SZ(graph[v]) == SZ(graph[u])))){
				for(int k=0; k < SZ(graph[v]); k++){
					int w = graph[v][k];

					if(w>u){																
						sign_sum += edge_sign[{v,w}];
						n_w.insert(w);
						
						if(sign_sum == 0) {count_wedge_with_signs_0[w] +=1;}
                        			else if(sign_sum == 1) {count_wedge_with_signs_1[w] +=1;}
                        			else if(sign_sum == 2) {count_wedge_with_signs_2[w] +=1;}	
					}
				}
			//	}
		}

		for(auto i : n_w)
		{
			            int two_zeros = count_wedge_with_signs_0[i];//- -            
            			    int two_ones = count_wedge_with_signs_2[i];// + +            
                                    int one_zero = count_wedge_with_signs_1[i];// - +

                                                                            
                                      if((two_zeros + two_ones)>2)
                                      res += ((two_zeros + two_ones) * ((two_zeros + two_ones) - 1) * ((two_zeros + two_ones) - 2))/6;
                                        if((one_zero)>2)
                                      	res += ((one_zero) * ((one_zero) - 1) * ((one_zero) - 2)) /6;                         
      
                                    
		}

	}

	return res;
}*/

  /*ll all_bal_bat_counting(vector<vector<int>> &graph) {
    ll result_count = 0;
    vector<int> count_wedge_0(n_vertices, 0);
    vector<int> count_wedge_1(n_vertices, 0);
    //vector<int> count_wedge_2(n_vertices, 0);
    vector<int> aux_array_two_neighboorhood(n_vertices);
    
    for (int i = 0; i < vertices_in_left.size(); i++) {
        int idx = 0;
        for (int j = 0; j < SZ(graph[i]); j++) {
            int v = graph[i][j];
            int sign_sum = edge_sign[{i, v}];
            
            for (int k = 0; k < SZ(graph[v]); k++) {
                int w = graph[v][k];
                
                if (w < i) {
                    sign_sum += edge_sign[{v, w}];
                    
                    if (sign_sum %2 == 0) {
                        count_wedge_0[w]++;
                    } else if (sign_sum %2 == 1) {
                        count_wedge_1[w]++;
                    } 
                    
                    if (count_wedge_0[w] == 1 || count_wedge_1[w] == 1) {
                        aux_array_two_neighboorhood[idx++] = w;
                    }
                } else {
                    break;
                }
            }
        }

        for (int j = 0; j < idx; j++) {
            int w = aux_array_two_neighboorhood[j];
            int total_similar = count_wedge_0[w];
            int total_non_similar = count_wedge_1[w];
            
            if (total_similar > 2) {
                result_count += (total_similar * (total_similar - 1) * (total_similar - 2)) / 6;
            }
            
            if (total_non_similar > 2) {
                result_count += (total_non_similar * (total_non_similar - 1) * (total_non_similar - 2)) / 6;
            }
            
            count_wedge_0[w] = 0;
            count_wedge_1[w] = 0;
            
        }
    }
    return result_count;
}*/


ll all_bal_bat_counting(vector < vector <int> > &graph) {
	ll res=0;
	ll balanced_bf_count=0;

	for(int u=0; u < vertices_in_left.size(); u++){
                unordered_map<int, int> count_wedge_with_signs_0; 
		unordered_map<int, int> count_wedge_with_signs_1; 
		 
		set<int> n_w;
		for(int j = 0; j < SZ(graph[u]); j++)
		{
			int v = graph[u][j];
			
			//if(SZ(graph[v]) < SZ(graph[u]) || ((SZ(graph[v]) == SZ(graph[u])) &&(v<u))){
				int sign_sum = edge_sign[{u,v}]; 
				int temp = sign_sum; 

				for(int k=0; k < SZ(graph[v]); k++){
					int w = graph[v][k];
					sign_sum = temp; 
					 
					if(SZ(graph[w]) < SZ(graph[u]) ||  ((SZ(graph[w]) == SZ(graph[u])) &&(w<u))) {
						sign_sum += edge_sign[{v,w}]; 
						n_w.insert(w); 
                        			if(sign_sum % 2 == 0) {count_wedge_with_signs_0[w] +=1;}
                        			else if(sign_sum % 2 == 1) {count_wedge_with_signs_1[w] +=1;}
                        			

					}
				}
			//}
		}
	
	
        for(auto i : n_w)
		{
            int two_zeros = count_wedge_with_signs_0[i];//- -
            int one_zero = count_wedge_with_signs_1[i];// - +
          
			if(two_zeros > 1)
			{
			   balanced_bf_count += (((two_zeros) * (two_zeros -1 ) * (two_zeros - 2)) /6);
			}
			
			if(one_zero > 1)
			{
			   balanced_bf_count += (((one_zero) *  (one_zero -1 )  * (one_zero - 2)) /6);
			}

			
		}

	}
	return balanced_bf_count;
}
void exact_algorithm_time_tracker(){

	double beg_clock, end_clock, elapsed_time;
	
	beg_clock = clock();    	
    	ll total_bat_count1 = all_bal_bat_counting(adj); 
	end_clock = clock();
    elapsed_time = (end_clock - beg_clock) / CLOCKS_PER_SEC;
     cout << "total time " << elapsed_time << " seconds." << endl;	
		
	printf("Total number of bats : %lld \n", total_bat_count1);
	
	
		
}


int main(int argc, char *argv[])
{
	std::ios::sync_with_stdio(false);
	read_the_graph(argv[1]);
	exact_algorithm_time_tracker();
}


