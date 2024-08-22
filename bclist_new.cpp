

#include "bclist_new.h"

VertexDegree::VertexDegree(): vertex(-1), degree(-1){}

VertexDegree::VertexDegree(int v, int d): vertex(v), degree(d){}

VertexDegree::~VertexDegree(){}



lint Tools::choose(lint n, lint k){
    if(combination_cache[n][k] > 0){
        return combination_cache[n][k];
    }
    else{
        if (k > n) {
            return 0;
        }
        lint backup_n = n;
        lint backup_k = k;
        lint r = 1;
        for (lint d = 1; d <= k; ++d) {
            r *= n--;
            r /= d;
        }
        if(n <= MAX_N && k <= MAX_K){
            combination_cache[backup_n][backup_k] = r;
        }
        return r;
    }
}

vector<int> Tools::toVector(const string &str){
    vector<int> ans;
    string c = "/";
    string::size_type pos1, pos2;
    pos2 = str.find(c);
    pos1 = 0;
    while(std::string::npos != pos2)
    {
        ans.emplace_back(stoi(str.substr(pos1, pos2-pos1)));
    
        pos1 = pos2 + c.size();
        pos2 = str.find(c, pos1);
    }
    if(pos1 != str.length())
        ans.emplace_back(stoi(str.substr(pos1)));
    
    return ans;
}

// below are special bigraphs

SpecialBigraph::SpecialBigraph(string in_path, int p_value, int q_value, int priority_strategy, int use_cost_model_flag){
    clear_everything();
    data_file_path = in_path;
    p = p_value;
    q = q_value;
    priority = priority_strategy;
    use_cost_model = use_cost_model_flag;
}

void SpecialBigraph::clear_everything(){
    largest_index_in_partition[0] = largest_index_in_partition[1] = 0;
    num_vertices = 0;
    num_edges = 0;
    result_count = 0;
    bcliques.clear();
    edges.clear();
    vertices[0].clear(); vertices[1].clear();
    vertices_in_left.clear();
    vertices_in_right.clear();
    adj_vec.clear();
}



void SpecialBigraph::add_vertex(int A, int side) {
    if(vertices[side].find(A) == vertices[side].end()) {
        if(side == 0)
            vertices_in_left.emplace_back(A);
        else
            vertices_in_right.emplace_back(A);
        vertices[side][A] = 1;
    }
}

void SpecialBigraph::add_edge(int &A, int &B, int &sign) {
    add_vertex(A, 0);
    add_vertex(B, 1);
    if (edges.find({A, B, sign}) == edges.end()) {
        edges.insert({A, B, sign});
        num_edges++;
    }
}

void SpecialBigraph::get_index(int &A, int side) {
    if(vertices[side].find(A) == vertices[side].end()) {
        vertices[side][A] = largest_index_in_partition[side]++;
    }
    A = vertices[side][A];
}

bool SpecialBigraph::all_num(string &s) {
    for (int i = 0; i < SZ(s); i++) 
        if ((s[i] >= '0' && s [i] <= '9') == false) return false;
    return true;
}

void SpecialBigraph::read_graph() {
    freopen(data_file_path.c_str(), "r", stdin);
    string s;
    cin.clear();
    while(getline(cin, s)){
        stringstream ss;
        ss << s;
        vector<string> vec_str;
        for(string z; ss >> z; vec_str.push_back(z));
        if(SZ(vec_str) >= 2){
            bool is_all_number = true;
            for(int i = 0; i < min(2, SZ(vec_str)); i++)
                is_all_number &= all_num(vec_str[i]);
            if(is_all_number){
                int A, B;
                int sign;
                ss.clear(); ss << vec_str[0]; ss >> A;
                ss.clear(); ss << vec_str[1]; ss >> B;
                ss.clear(); ss << vec_str[2]; ss >> sign;
                add_edge(A, B, sign);
            }
        }
    }

    vertices[0].clear();
    vertices[1].clear();
    
    n_vertices[0] = vertices_in_left.size(); 
    n_vertices[1] = vertices_in_right.size();
    
    num_vertices = n_vertices[0] + n_vertices[1];
    
    num_edges = edges.size();
    deg = new int[num_vertices](); // Initialize degrees to 0
    adj_vec.resize(num_vertices, vector<int>());

    largest_index_in_partition[0] = 0;
    largest_index_in_partition[1] = n_vertices[0];    
    adj_vec.resize(num_vertices, vector<int>());
    for (auto edge : edges) {
        int A = edge[0];  
        int B = edge[1];
        int sign = edge[2];
       
        get_index(A, 0);
        get_index(B, 1);
        
        adj_vec[A].push_back(B);
        adj_vec[B].push_back(A);
               
               
          // Increment the degree for both vertices
        deg[A]++;
        deg[B]++;     
                 
        edge_sign[{A, B}] = sign;
        edge_sign[{B, A}] = sign;    
        
        //cout << "Edge (" << A << ", " << B << ") has sign: " << sign << endl;
        
    }
    

/*    #ifdef IS_DEBUGGING
    print_adj();
    print_edges();
    #endif*/
    
    cout << "end of read graph" << endl;
    cout << "#vertices = " << num_vertices << ", #left_vertices = " << n_vertices[0] << ", #right_vertices = " << n_vertices[1] << endl;
    cout << "#edges = " << num_edges << endl;
}


void SpecialBigraph::trim_graph_by_core(){

    int start_idx, end_idx;
    start_idx = end_idx = 0;
    
    int *to_remove_vertices = new int[num_vertices];
    bool *removed = new bool[num_vertices];
    
    
    for(int i = 0; i < num_vertices; i++){
        to_remove_vertices[i] = 0;
        removed[i] = false;
   
    }
    
    // step 1: collect the initial vertices with degree less than q in left and p in right
    for(int i = 0; i < num_vertices; i++){
        if((i < n_vertices[0] && deg[i] < q) || (i >= n_vertices[0] && deg[i] < p)){
            to_remove_vertices[end_idx++] = i;
            removed[i] = true;
            
        }
    }    

    // step 2: recursively remove all vertices with degree less than q in left and p in right, i.e., (q,p)-core
    while(start_idx != end_idx){
        int vertex = to_remove_vertices[start_idx++];
        
        //cout << "remove : " << vertex << endl;
        for(int i = 0; i < adj_vec[vertex].size(); i++){
            int other = adj_vec[vertex][i];
            if(!removed[other]){
                deg[other]--;
                
                if((other < n_vertices[0] && deg[other] < q) || (other >= n_vertices[0] && deg[other] < p)){
                    to_remove_vertices[end_idx++] = other;
                    removed[other] = true;
                }
            }
        }
        adj_vec[vertex].clear();
        deg[vertex] = 0;
    }
    
    

    for(int i = 0; i < num_vertices; i++){
        if(!removed[i]){
            vector<int> new_vec(deg[i]);
            int idx = 0;
            for(int j = 0; j < adj_vec[i].size(); j++){
                if(!removed[adj_vec[i][j]]){
                    new_vec[idx++] = adj_vec[i][j];
                }
            }
            adj_vec[i].clear();
            adj_vec[i] = new_vec;
            
  
        }
    }
    delete[] to_remove_vertices;
    delete[] removed;
    
    #ifdef IS_DEBUGGING

    print_adj();
    #endif
    
    cout << "after trimming by core" << endl;
    
    reformat_graph();
    
    cout << "end of trim graph by core" << endl;
}

void SpecialBigraph::reformat_graph(){
    // re-format the bipartite graph since we might have removed some vertices
    edges.clear();
    for(int i = 0; i < n_vertices[0]; i++){
        for(int j = 0; j < deg[i]; j++){
            //edges.insert(make_pair(i, adj_vec[i][j]));
           edges.insert({i, adj_vec[i][j], edge_sign[{i, adj_vec[i][j]}]});


        }
    }
    num_edges = edges.size();
    
    int removed_count[2];
    removed_count[0] = removed_count[1] = 0;
    int side;
    for(int i = 0; i < num_vertices; i++){
        if(deg[i] == 0){
            side = i < n_vertices[0] ? 0 : 1;
            removed_count[side]++;
        }
    }
    
    n_vertices[0] = n_vertices[0] - removed_count[0];
    n_vertices[1] = n_vertices[1] - removed_count[1];
    
    num_vertices = n_vertices[0] + n_vertices[1];
    
    list_of_edges.clear();
    for(int i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    adj_vec.clear();
    
    //cout << num_vertices << endl;
    
    vertices[0].clear();
    vertices[1].clear();
    delete[] deg;
    deg = new int[num_vertices];
    memset(deg, 0, num_vertices * sizeof(int));
    adj_vec.resize(num_vertices, vector<int>());
    
    largest_index_in_partition[0] = 0;
    largest_index_in_partition[1] = n_vertices[0];
    

for (auto edge : edges) {
        int A = edge[0];  
        int B = edge[1];
        int sign = edge[2];
       
        get_index(A, 0);
        get_index(B, 1);
        
        adj_vec[A].push_back(B);
        adj_vec[B].push_back(A);
               
               
          // Increment the degree for both vertices
        deg[A]++;
        deg[B]++;     
                 
        edge_sign[{A, B}] = sign;
        edge_sign[{B, A}] = sign;    
        
        cout << "Edge (" << A << ", " << B << ") has sign: " << sign << endl;
        
    }

    
    for(int i = 0; i < num_vertices; i++){
        sort(adj_vec[i].begin(), adj_vec[i].end());
    }
    
    #ifdef IS_DEBUGGING
    cout<<"=============print adj=========="<<endl;
    print_adj();  
    cout<<"=============print edges=========="<<endl;
    print_edges();    
    #endif
    
    cout << "end of reformat graph" << endl;
    
    cout << "#vertices = " << num_vertices << ", #left_vertices=" << n_vertices[0] << ", #right_vertices=" << n_vertices[1] << endl;
    cout << "#edges = " << edges.size() << endl;
   
}


void SpecialBigraph::sort_vertices(int strategy){
    switch (strategy) {

        case 1: {// degree vertex order
            cout << "degree vertex order is applied..." << endl;
            
            //cout << "num_vertices = " << num_vertices << endl;
            
            vector<VertexDegree> verdegs(num_vertices);
            for(int i = 0; i < num_vertices; i++){
                VertexDegree vd(i, deg[i]);
                verdegs[i] = vd;
            }
            sort(verdegs.begin(), verdegs.end());
            
            all_vertices.resize(num_vertices);
            for(int i = 0; i < num_vertices; i++){
                all_vertices[i] = verdegs[i].vertex;
            }
            
            break;
        }
        default:{
            cout << "Please select the following vertex ordering: 0-random, 1-degree, 2-core" << endl;
            return;
        }
    }
    vector<int> left_vertices_tmp(n_vertices[0]);
    vector<int> right_vertices_tmp(n_vertices[1]);
    int left_id = 0;
    int right_id = 0;
    for(int i = 0; i < all_vertices.size(); i++){
        if(all_vertices[i] < n_vertices[0]){
            left_vertices_tmp[left_id++] = all_vertices[i];
        }
        else{
            right_vertices_tmp[right_id++] = all_vertices[i];
        }
    }
    
    all_vertices.clear();
    if(anchor_left){
        all_vertices = left_vertices_tmp;
        all_vertices.insert(all_vertices.end(), right_vertices_tmp.begin(), right_vertices_tmp.end());
    }
    else{
        all_vertices = right_vertices_tmp;
        all_vertices.insert(all_vertices.end(), left_vertices_tmp.begin(), left_vertices_tmp.end());
    }
}


void SpecialBigraph::collect_two_hop_adj(){
    int array_size = n_vertices[0];
    int array_sizen = n_vertices[1];
    cout<<"array_size :"<<array_size<<endl;
    cout<<"array_sizen :"<<array_sizen<<endl;
    cout<<"================================================================="<<endl;
    two_hop_adj_maxsize = new int[array_size];
    two_hop_adj_vec = new int*[array_size];
    for(int i = 0; i < array_size; i++){
        two_hop_adj_maxsize[i] = NTWOHOPS;
        two_hop_adj_vec[i] = new int[two_hop_adj_maxsize[i]];
    }
    
    two_hop_adj_size = new int[array_size]();
    int *common_neig_map = new int[array_size]();
    int *aux_array_two_neig = new int[array_size]();
    
    for(int i = 0; i < array_size; i++){
        int idx = 0;       
        for(int j = 0; j < adj_vec[i].size(); j++){
            int v = adj_vec[i][j];
            for(int k = 0; k < adj_vec[v].size(); k++){
                int w = adj_vec[v][k];
                if(w < i){
                cout<<"U "<<i<<endl;
                cout<<"V "<<v<<endl;
                cout<<"W "<<w<<endl;
                    common_neig_map[w]++;
                    if(common_neig_map[w] == 1){
                        aux_array_two_neig[idx++] = w;
                    }
                }
                else{
                    break;
                }
            }
        }
        for(int j = 0; j < idx; j++){
            int w = aux_array_two_neig[j];
            if(common_neig_map[w] >= q){            
                if(two_hop_adj_size[w] >= two_hop_adj_maxsize[w]){                  
         
                    two_hop_adj_maxsize[w] *= 2;
                    int *temp_array = new int[two_hop_adj_maxsize[w]];
                    memcpy(temp_array, two_hop_adj_vec[w], sizeof(int)*two_hop_adj_size[w]);
                    delete[] two_hop_adj_vec[w];
                    two_hop_adj_vec[w] = temp_array;
                }
                two_hop_adj_vec[w][two_hop_adj_size[w]++] = i;
            }
            common_neig_map[w] = 0;            
        }                        
    }

    delete[] common_neig_map;
    delete[] aux_array_two_neig;

}

bool SpecialBigraph::prepare_graph(){
    #ifdef USE_CORE_REDUCTION
    trim_graph_by_core();
    #endif
    
    if(n_vertices[0] == 0 || n_vertices[1] == 0){
        cout << "No results because the graph is pruned by core" << endl;
        return false;
    }
    
    if(estimate_cost(0) < estimate_cost(1)){
        anchor_left = true;
    }
    else{
        anchor_left = false;
    }
    
    if(use_cost_model == 0){
        if(anchor_left == true)
            anchor_left = false;
        else
            anchor_left = true;
    }
    
    cout << "anchor_left = " << anchor_left << endl;  

    
    sort_vertices(priority);
    
    #ifdef IS_DEBUGGING
    print(all_vertices);
    #endif
    cout << "finish sorting vertices" << endl;
  
    int id = 0;
    for(vector<int>::iterator it = all_vertices.begin(); it != all_vertices.end(); it++){
        vertexids[*it] = id++;
    }
	
    
    // sort the graph by the rank of vertices
    for(int i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    

    cout << "========================sort the graph by the rank of vertices=======================" << endl;

    #ifdef IS_DEBUGGING
    print(all_vertices);
    #endif
    
    #ifdef IS_DEBUGGING
    cout<<"=============print adj=========="<<endl;
    print_adj();  
    cout<<"=============print edges=========="<<endl;
    print_edges();    
    #endif    
    
     cout<<"=====================  adj_vec. clear() ============================================"<<endl;
    adj_vec.clear();
        #ifdef IS_DEBUGGING
    cout<<"=============print adj=========="<<endl;
    print_adj();  
    cout<<"=============print edges=========="<<endl;
    print_edges();    
    #endif   
    delete[] deg;
    deg = new int[num_vertices]();
    adj_vec.resize(num_vertices, vector<int>());

   for (auto edge : edges) {
        
        int A = edge[0];  
        int B = edge[1];
        int sign = edge[2];
       
        get_index(A, 0);
        get_index(B, 1);
        
        adj_vec[A].push_back(B);
        adj_vec[B].push_back(A);                                       
        deg[A]++;
        deg[B]++;     
                 
        edge_sign[{A, B}] = sign;
        edge_sign[{B, A}] = sign;    
        
    }		

       if(!anchor_left){
       //cout<<"anchor_left "<<anchor_left<<endl;
       
        int tmp_num = n_vertices[0];
        n_vertices[0] = n_vertices[1];
        n_vertices[1] = tmp_num;
        
        int tmp_value = p;
        p = q;
        q = tmp_value;
    }
    cout<<"p " << p  <<endl;
    cout<<"q " << q <<endl;
    for(int i = 0; i < num_vertices; i++){
        sort(adj_vec[i].begin(), adj_vec[i].end());
    }
    

    collect_two_hop_adj();
   for(int i = 0; i < n_vertices[0]; i++){
        sort(two_hop_adj_vec[i], two_hop_adj_vec[i] + two_hop_adj_size[i]);
    }
    
    #ifdef IS_DEBUGGING
    cout << "adj after preparing" << endl;
    print_adj();
    print_two_hop_adj();
    cout<<"----------------------------------------------------"<<endl;
    #endif
    return true;
}

double SpecialBigraph::estimate_cost(int side){
    srand (time(NULL));
    
    int num_rounds = ceil(num_vertices * SAMPLE_RATE);
    lint total_two_hop_deg = 0;
    lint max_two_hop_deg = 0;
    
    int *common_neig_map = new int[num_vertices]();
    int *aux_array_two_neig = new int[num_vertices]();
    int offset = side == 0 ? 0 : n_vertices[0];
    int common_neig_threshold = side == 0 ? q : p;
  
    for(int r = 0; r < num_rounds; r++){
        lint estimated_two_hop_deg = 0;
        int u = rand() % n_vertices[side] + offset;
        cout<<"u "<<u<<endl;
        int idx = 0;
        for(int j = 0; j < adj_vec[u].size(); j++){
            int v = adj_vec[u][j];
            for(int k = 0; k < adj_vec[v].size(); k++){
                int w = adj_vec[v][k];
                if(w < u){
                    common_neig_map[w]++;
                    if(common_neig_map[w] == 1){
                        aux_array_two_neig[idx++] = w;
                    }
                }
                else{
                    break;
                }
            }
        }		 
        for(int j = 0; j < idx; j++){
            int w = aux_array_two_neig[j]; // 0
            if(common_neig_map[w] >= common_neig_threshold){
                estimated_two_hop_deg += 1;
            }
            common_neig_map[w] = 0;
        }
        
        max_two_hop_deg = max_two_hop_deg < estimated_two_hop_deg ? estimated_two_hop_deg : max_two_hop_deg;       
        total_two_hop_deg = total_two_hop_deg + estimated_two_hop_deg * n_vertices[side];
    }
    
    total_two_hop_deg = total_two_hop_deg / num_rounds;
    cout << "estimated_total_two_hop_deg[" << side << "]=" << total_two_hop_deg << ", estimated_avg_two_hop_deg[" << side << "]=" << total_two_hop_deg / n_vertices[side] << ", max_two_hop_deg[" << side << "]=" << max_two_hop_deg << endl;
    
    lint avg_two_hop_deg = max(2, total_two_hop_deg / n_vertices[side]);  // we let avg_two_hop_deg be at least 2 to deal with corner case
    lint pq_value = side == 0? p:q;
    double totalCost = total_two_hop_deg * pow(avg_two_hop_deg, pq_value-2);
   cout << "totalCost = " << totalCost << endl;
    
    delete[] common_neig_map;
    delete[] aux_array_two_neig;
    return totalCost;
}
 
 
 
 void SpecialBigraph::listing_cliques(){
    int start_idx, max_d, max_two_hop_d, tmp_ns, e, two_hop_e;
    int *tmp_sub, *tmp_two_hop_d;
    
      n = n_vertices[0];
    
    d = new int[n];
    tmp_two_hop_d = new int[n];
    
    e = two_hop_e = 0;
    for(int i = 0; i < n; i++){
        d[i] = adj_vec[i].size();
        e += d[i];
    }
    
    for(int i = 0; i < n; i++){      
        tmp_two_hop_d[i] = two_hop_adj_size[i];
        two_hop_e += tmp_two_hop_d[i];
    }

    ns = new int[p+1]; // 3
    cd = new int[n+1]; // 3
    adj = new int[e]; // 6
    two_hop_d = new int*[p+1]; //3
    two_hop_cd = new int[n+1]; // 3
    two_hop_adj = new int[two_hop_e]; // 1
    lab = new int[n]; 
    sub = new int*[p+1];
    op_vertices = new int*[p+1];
    op_size = new int[p+1];
    result_count = 0;
	
    
    tmp_ns = 0;
    cd[0] = 0;
    two_hop_cd[0] = 0;
    max_d = 0;
    max_two_hop_d = 0;
    tmp_sub = new int[n]; 
    
    for(int i = 1; i < n+1; i++){
        cd[i] = cd[i-1] + d[i-1];
        max_d = (max_d > d[i-1])?max_d:d[i-1];
        two_hop_cd[i] = two_hop_cd[i-1] + tmp_two_hop_d[i-1];
        max_two_hop_d = (max_two_hop_d > tmp_two_hop_d[i-1])?max_two_hop_d:tmp_two_hop_d[i-1];
        tmp_sub[tmp_ns++] = i-1;
        lab[i-1] = p;
        
        
        // Print the values for each iteration
        std::cout << "Iteration " << i << ":\n";
        std::cout << "cd[" << i << "] = " << cd[i] << "\n";
        std::cout << "max_d = " << max_d << "\n";
        std::cout << "two_hop_cd[" << i << "] = " << two_hop_cd[i] << "\n";
        std::cout << "max_two_hop_d = " << max_two_hop_d << "\n";
        std::cout << "tmp_sub[" << tmp_ns-1 << "] = " << tmp_sub[tmp_ns-1] << "\n";
        std::cout << "lab[" << i-1 << "] = " << lab[i-1] << "\n";
         std::cout << "tmp_ns " << tmp_ns << ":\n";
        std::cout << "----------------------------\n";

    }
 
	     for(int i = 0; i < n; i++){
		for(int j = 0; j < adj_vec[i].size(); j++){
		    adj[cd[i] + j] = adj_vec[i][j];
		}

    	     }
    	     
	    for(int i = 0; i < n; i++){	   
		for(int j = 0; j < two_hop_adj_size[i]; j++){
		    two_hop_adj[two_hop_cd[i] + j] = two_hop_adj_vec[i][j];
		}
	    }

     
for(int i = 1; i < p; i++) {
    two_hop_d[i] = new int[n];
    sub[i] = new int[max_two_hop_d];
    op_vertices[i] = new int[max_d];

    // Print the indices and the sizes of the allocated arrays
    std::cout << "Index i: " << i << std::endl;
    std::cout << "Allocated two_hop_d[" << i << "] with size " << n << std::endl;
    std::cout << "Allocated sub[" << i << "] with size " << max_two_hop_d << std::endl;
    std::cout << "Allocated op_vertices[" << i << "] with size " << max_d << std::endl;
    std::cout << std::endl;  // Print an empty line for better readability
}


      		ns[p] = tmp_ns;

    		two_hop_d[p] = tmp_two_hop_d;
    		
    		sub[p] = tmp_sub;
    		
	
		op_vertices[p] = new int[max_d];

           cout<<"-----------------------------------------------"<<endl;
       
	    cout << "construct graph done" << endl;
	    cout << "Calling pqclique with p = " << p << endl;  
	    clique_vertices_in_left.resize(p);
	     pqclique(p);
    
}
 
void SpecialBigraph::pqclique(int l){
    //cout << "ns[" << l << "]=" << ns[l] << endl;
    
    int a,i,j,k,end,u,v,w,sign_sum ,sign_index=0,sign_index1=0;
    std::vector<int> sign_vec;
    std::vector<int> sign_vec1;
    sign_vec.resize(3);
    sign_vec1.resize(3); // Vector to store the signs
    if(l == 2){
    cout<<"-----------------------------------------------"<<endl;
        cout<<"ns[2] "<<ns[2]<<endl; //2
        
        // Initialize the sign index counter
        sign_index = 0;
        for(i = 0; i < ns[2]; i++){
            u = sub[2][i];
            cout<<"u " << u<<endl; // 0
            clique_vertices_in_left[p-2] = u;
            
            op_size[2] = 0;
            end = cd[u] + d[u]; // 3
            if(p <= 2){
                for(j = cd[u]; j < end; j++){
// op_vertices[2][op_size[2]++] = adj[j];
			op_vertices[2][op_size[2]++] = adj[j];
			sign_sum = edge_sign[{u,adj[j]}];
			// Store the sign in the vector sequentially
                         sign_vec[sign_index++] = sign_sum;
			// Print the sign value and its index in the vector
                    std::cout << "sign_vec[" << sign_index - 1 << "] = " << sign_sum << std::endl;
                }
            }
            else{
                j = 0;
                k = cd[u];
                while(j < op_size[3] && k < cd[u] + d[u]){
                    if(op_vertices[3][j] == adj[k]){
                        op_vertices[2][op_size[2]++] = adj[k];
                        j++;
                        k++;
                    }
                    else if(op_vertices[3][j] < adj[k]){
                        j++;
                    }
                    else{
                        k++;
                    }
                }
            }
            
            if(op_size[l] < q){
                continue;
            }
            
            end = two_hop_cd[u] + two_hop_d[2][u];
            cout<<"end " <<end<<endl;
            for(a = two_hop_cd[u]; a < end; a++){
                v = two_hop_adj[a];
                cout<<"v " <<v<<endl;                
                clique_vertices_in_left[p-1] = v;
                
                op_size[1] = 0;
                
                j = 0;
                k = cd[v];
                while(j < op_size[2] && k < cd[v] + d[v]){
                    if(op_vertices[2][j] == adj[k]){                     
                        op_vertices[1][op_size[1]++] = adj[k];
                        sign_sum = edge_sign[{v,adj[k]}];		
                        sign_vec1[sign_index1++] = sign_sum;		
                    	std::cout << "sign_vec1[" << sign_index1 - 1 << "] = " << sign_sum << std::endl;
                        j++;
                        k++;
                    }
                    else if(op_vertices[2][j] < adj[k]){
                        j++;
                    }
                    else{
                        k++;
                    }
                }
                if(op_size[1] < q){
                    continue;
                }


int equal_count = 0;
int not_equal_count = 0;
int two_three_bal = 0;

// Compare sign_vec and sign_vec1 and count equal and non-equal pairs
for(int i = 0; i < sign_vec.size(); i++) {
    if(sign_vec[i] == sign_vec1[i]) {
        equal_count++;
    } else {
        not_equal_count++;
    }
}

// Perform "choose 3" operation for equal and non-equal counts
int equal_combination = (equal_count >= 3) ? (equal_count * (equal_count - 1) * (equal_count - 2)) / 6 : 0;
int not_equal_combination = (not_equal_count >= 3) ? (not_equal_count * (not_equal_count - 1) * (not_equal_count - 2)) / 6 : 0;

// Output the results
std::cout << "C(equal_count, 3) = " << equal_combination << std::endl;
std::cout << "C(not_equal_count, 3) = " << not_equal_combination << std::endl;

two_three_bal = equal_combination + not_equal_combination;
std::cout << "two_three_bal = " << two_three_bal << std::endl;

                
            }
            
        }
        return;
    }
    for(i = 0; i < ns[l]; i++){
        u = sub[l][i];
        //clique_vertices_in_left[l-1] = u;
        clique_vertices_in_left[p-l] = u;
        //cout << "u[" << l << "]=" << u << endl;
        // to handle opposite vertices;
        op_size[l] = 0;
        if(l == p){
            // op_vertices[p] is simply all adj vertices of u
            end = cd[u] + d[u];
            for(j = cd[u]; j < end; j++){
                op_vertices[l][op_size[l]++] = adj[j];
            }
        }
        else{
            j = 0;
            k = cd[u];
            while(j < op_size[l+1] && k < cd[u] + d[u]){
                if(op_vertices[l+1][j] == adj[k]){
                    op_vertices[l][op_size[l]++] = adj[k];
                    j++;
                    k++;
                }
                else if(op_vertices[l+1][j] < adj[k]){
                    j++;
                }
                else{
                    k++;
                }
            }
        }
        
        if(op_size[l] < q){
            continue;
        }
        
        ns[l-1] = 0;
        end = two_hop_cd[u] + two_hop_d[l][u];
        //cout << "two_hop_cd[" << u << "] = " << two_hop_cd[u] << endl;
        //cout << "two_hop_d[" << l << "][" << u << "] = " << two_hop_d[l][u] << endl;
        for(j = two_hop_cd[u]; j < end; j++){
            v = two_hop_adj[j];     //cout << "v=" << v << endl;
            if(lab[v] == l){
                lab[v] = l-1;
                sub[l-1][ns[l-1]++] = v;
                two_hop_d[l-1][v] = 0;
            }
        }
        for(j = 0; j < ns[l-1]; j++){
            v = sub[l-1][j];
            end = two_hop_cd[v] + two_hop_d[l][v];
            for(k = two_hop_cd[v]; k < end; k++){
                w = two_hop_adj[k];
                if(lab[w] == l-1){
                    two_hop_d[l-1][v]++;
                }
                else{
                    two_hop_adj[k--] = two_hop_adj[--end];
                    two_hop_adj[end] = w;
                }
            }
        }
        
        pqclique(l-1);
        
        for(j = 0; j < ns[l-1]; j++){
            v = sub[l-1][j];
            lab[v] = l;
        }
    }
}
  

void SpecialBigraph::count_butterflies() {

    if (estimate_cost(0) < estimate_cost(1)) {
        anchor_left = true;
    } else {
        anchor_left = false;
    }
    
    cout << "use_cost_model " << use_cost_model << endl;

    if (use_cost_model == 0) {
        anchor_left = !anchor_left;
    }
    
    cout << "anchor_left = " << anchor_left << endl;

    int start_idx = anchor_left ? n_vertices[0] : 0;
    int end_idx = anchor_left ? num_vertices : n_vertices[0];

    for (int i = start_idx; i < end_idx; i++) {
        sort(adj_vec[i].begin(), adj_vec[i].end());
    }
      
    start_idx = anchor_left ? 0 : n_vertices[0];
    end_idx = anchor_left ? n_vertices[0] : num_vertices;
    
    
    cout<<"-------------------------graph after reading----------------------- "<<endl;
    
    #ifdef IS_DEBUGGING
    print_adj();
    print_edges();
    #endif
    cout<<"-------------------------------------------------------------------- "<<endl;
ll balanced_bf_count=0;
	for(int u = start_idx; u < end_idx; u++){
		set<int> n_w;
		unordered_map<int,vector<int>> wedges;
		for(int j = 0; j < adj_vec[u].size(); j++)
		{
			int v = adj_vec[u][j];
				for(int k=0; k < adj_vec[v].size(); k++){
					int w = adj_vec[v][k];
					 
					if(w < u){									
						wedges[w].push_back(v);
						n_w.insert(w); 
					}
				}			
		}
		
		for(auto i : n_w){
			vector<int> vs = wedges[i];
			ll u1 = u;
			ll w1 = i;
			for(int j=0; j < vs.size()-1; j++){
				ll v1 = vs[j];

				for(int k=j+1; k < vs.size();k++){
					ll x1 = vs[k];
					int sign_sum1 = edge_sign[{u1,x1}] + edge_sign[{w1,x1}] +edge_sign[{u1,v1}] + edge_sign[{w1,v1}];
					if(sign_sum1 == 0 || sign_sum1 == 2 || sign_sum1 == 4) balanced_bf_count++;
				}
				

			}
		}
        
	}

	cout<<" balanced_bf_count: "<< balanced_bf_count<<endl; 
}

    
SpecialBigraph::~SpecialBigraph(){
    delete[] deg;
    
    vertices[0].clear();
    vertices[1].clear();
    vertices_in_left.clear();
    vertices_in_right.clear();
    for(int i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    adj_vec.clear();
    for(int i = 0; i < n_vertices[0]; i++){
        delete[] two_hop_adj_vec[i];
    }
    delete[] two_hop_adj_vec;
    delete[] two_hop_adj_size;
    delete[] two_hop_adj_maxsize;
    edges.clear();
    list_of_edges.clear();
    
    delete[] ns;
    delete[] d;
    delete[] cd;
    delete[] adj;
    for(int i = 1; i <= p; i++){
        delete[] two_hop_d[i];
        delete[] sub[i];
        delete[] op_vertices[i];
    }
    delete[] two_hop_d;
    delete[] two_hop_cd;
    delete[] two_hop_adj;
    delete[] lab;
    delete[] sub;
    delete[] op_vertices;
    delete[] op_size;
}


// the following are tool functions

void SpecialBigraph::print(vector<int> vec){
    for(vector<int>::iterator it = vec.begin(); it != vec.end(); it++){
        cout << *it << " ";
    }
    cout << endl;
}

void SpecialBigraph::print_adj(){
    for(int i = 0; i < num_vertices; i++){
        cout << "adj[" << i << "] : ";
        for(int j = 0; j < deg[i]; j++){
            cout << adj_vec[i][j] << " ";
        }
        cout << endl;
    }
}

void SpecialBigraph::print_adj(int i){
    for(int j = 0; j < deg[i]; j++){
        cout << adj_vec[i][j] << " ";
    }
    cout << endl;
}

void SpecialBigraph::print_two_hop_adj(){
    
    int maxDeg = 0;
    
    for(int i = 0; i < n_vertices[0]; i++){
        maxDeg = maxDeg >= two_hop_adj_size[i] ? maxDeg : two_hop_adj_size[i];
        
        cout << "two_hop_adj[" << i << "] : ";
        for(int j = 0; j < two_hop_adj_size[i]; j++){
            cout << two_hop_adj_vec[i][j] << " ";
        }
        cout << endl;
    }
    cout << "max two hop degree : " << maxDeg << endl;
}

/*void SpecialBigraph::print_edges(){
    cout << "#vertices : " << num_vertices <<endl;
    cout << "#edges : " << num_edges << endl;

}*/


void SpecialBigraph::print_edges() {
    cout << "#vertices : " << num_vertices <<endl;
    cout << "#edges : " << num_edges << endl;
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < deg[i]; j++) {
            int neighbor = adj_vec[i][j];
            // Since this is a bipartite graph, you may want to print edges only once
            if (i < neighbor) {
                cout << i << "  " << neighbor << endl;
            }
            

        }
    }
}


void SpecialBigraph::print_deg(){
    sort(deg, deg+num_vertices);
    for(int i = num_vertices - 1; i >= num_vertices - 1000; i--){
        cout << deg[i] << endl;
    }
}

void SpecialBigraph::print_bclique(ostream &os, vector<int> &left_vertices, vector<int> &right_vertices){
    for(int i = 0; i < left_vertices.size(); i++){
        #ifdef LARGE_BICLIQUES
        os << original_ids[left_vertices[i]] << " ";
        #else
        os << left_vertices[i] << " ";
        #endif
    }
    os << "| ";
    for(int i = 0; i < right_vertices.size(); i++){
        #ifdef LARGE_BICLIQUES
        os << original_ids[right_vertices[i]] << " ";
        #else
        os << right_vertices[i] << " ";
        #endif
    }
    os << endl;
}

void SpecialBigraph::print_results(){
    #ifdef COUNT_ONLY
    cout << "Total # results : " << result_count << endl;
    #endif
}   
    
void SpecialBigraph::print_map(unordered_map<int, int> ids){
    cout << "new mapping is as follows" << endl;
    for(unordered_map<int, int>::iterator it = ids.begin(); it != ids.end(); it++){
        cout << it->first << ", " << it->second << endl;
    }
}
    

int main(int argc, char *argv[]){
    double beg_clock, end_clock, elapsed_time;
    
    printf( "argc=%d\n", argc );
    
    SpecialBigraph *sbgraph;
    if(argc < 5){
        cout << "Too few arguments" << endl;
        return 0;
    }
    else if(argc == 5){
        for( int i = 0; i < argc; ++i )
            printf( "argv[%d]=%s\n", i, argv[i] );
        sbgraph = new SpecialBigraph(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    }
    else{
        for( int i = 0; i < argc; ++i )
            printf( "argv[%d]=%s\n", i, argv[i] );
        sbgraph = new SpecialBigraph(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
    }
    
    beg_clock = clock();
    sbgraph->read_graph();
    
    

    
    end_clock = clock();
    elapsed_time = (end_clock - beg_clock) / CLOCKS_PER_SEC;
    cout << "Read graph: using " << elapsed_time << " seconds." << endl;
    
    beg_clock = clock();
    #ifdef COUNT_ONLY
    if(atoi(argv[2]) == 2 && atoi(argv[3]) == 2){ //counting butterflies
        sbgraph->count_butterflies();
        end_clock = clock();
        elapsed_time = (end_clock - beg_clock) / CLOCKS_PER_SEC;
        cout << "Count butterflies: using " << elapsed_time << " seconds." << endl;
        sbgraph->print_results();
        return 0;
    }
   else if(!sbgraph->prepare_graph())
        return 0;
    #else
    if(!sbgraph->prepare_graph())
        return 0;
    #endif  
    end_clock = clock();
    elapsed_time = (end_clock - beg_clock) / CLOCKS_PER_SEC;
    cout << "Construct graph: using " << elapsed_time << " seconds." << endl;  
    
    beg_clock = clock();
    sbgraph->listing_cliques();
    end_clock = clock();
    elapsed_time = (end_clock - beg_clock) / CLOCKS_PER_SEC;
    cout << "Computing biclique: using " << elapsed_time << " seconds." << endl;
    
    
    return 0;
}
