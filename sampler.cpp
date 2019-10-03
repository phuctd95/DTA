#include "sampler.h"
#include <algorithm>

Sampler::Sampler(graph &gr):g(gr)
{
    //ctor
}

Sampler::~Sampler()
{
    //dtor
}

void Sampler::init_seed(int32_t p_init, vector <uint64_t> &s_init)
{
    rand_gen.init_seed(p_init,s_init);
}

void Sampler::get_random_pair(int32_t &s, int32_t &t)
{
    s = rand_gen.get_max(g.get_n_vertices()) + 1;
    while(1)
    {
        t = rand_gen.get_max(g.get_n_vertices()) + 1;
        if (s != t) break;
    }
}

int32_t Sampler::get_random_vertex()
{
    return rand_gen.get_max(g.get_n_vertices()) + 1;
}

void Sampler::init()
{
    shortest_path_from_s = vector <int32_t> (g.get_n_vertices()+1,-1);
    shortest_path_from_t = vector <int32_t> (g.get_n_vertices()+1,-1);
    n_paths_from_s = vector <double> (g.get_n_vertices()+1,0);
    n_paths_from_t = vector <double> (g.get_n_vertices()+1,0);
    sum_deg_s = 0;
    sum_deg_t = 0;
    path.clear();
    layer.clear();
}

void Sampler::add_s(int32_t u)
{
    queue_s.push(u);
    sum_deg_s += g.get_deg(u);
}

void Sampler::add_t(int32_t u)
{
    queue_t.push(u);
    sum_deg_t += g.get_deg_rev(u);
}

int32_t Sampler::get_front_s()
{
    int32_t u = queue_s.front();
    queue_s.pop();
    sum_deg_s -= g.get_deg(u);
    return u;
}

int32_t Sampler::get_front_t()
{
    int32_t u = queue_t.front();
    queue_t.pop();
    sum_deg_t -= g.get_deg_rev(u);
    return u;
}

void Sampler::add_from_s()
{
    int32_t p = shortest_path_from_s[queue_s.front()];
    int32_t u,v;
    while ((!queue_s.empty()) && (shortest_path_from_s[queue_s.front()] == p))
    {
        u = get_front_s();
        neigh = g.get_adjlist(u);
        for (int32_t i = 0; i < neigh.size(); ++i)
        {
            v = neigh[i];
            if (shortest_path_from_s[v] == -1)
            {
                shortest_path_from_s[v] = p + 1;
                add_s(v);
                if (shortest_path_from_t[v] != -1)
                {
                    stop = true;
                    layer.push_back(v);
                }
            }
            if (shortest_path_from_s[v] == p + 1)
                n_paths_from_s[v] += n_paths_from_s[u];
        }
        neigh.clear();
    }
}

void Sampler::add_from_t()
{
    int32_t p = shortest_path_from_t[queue_t.front()];
    int32_t u,v;
    while ((!queue_t.empty()) && (shortest_path_from_t[queue_t.front()] == p))
    {
        u = get_front_t();
        neigh = g.get_adjlist_rev(u);
        for (int32_t i = 0; i < neigh.size(); ++i)
        {
            v = neigh[i];
            if (shortest_path_from_t[v] == -1)
            {
                shortest_path_from_t[v] = p + 1;
                add_t(v);
                if (shortest_path_from_s[v] != -1)
                {
                    stop = true;
                    layer.push_back(v);
                }
            }
            if (shortest_path_from_t[v] == p + 1)
                n_paths_from_t[v] += n_paths_from_t[u];
        }
        neigh.clear();
    }
}

void Sampler::find_shortest_path(int32_t &s, int32_t &t)
{
    init();
    shortest_path_from_s[s] = 0;
    n_paths_from_s[s] = 1;
    shortest_path_from_t[t] = 0;
    n_paths_from_t[t] = 1;
    add_s(s);
    add_t(t);
    stop = false;
    while (!stop)
    {
        if ((queue_s.empty()) || (queue_t.empty())) break;
//            add_from_s(g);
        if (sum_deg_s <= sum_deg_t) add_from_s();
            else add_from_t();
    }
    while (!queue_s.empty()) queue_s.pop();
    while (!queue_t.empty()) queue_t.pop();
}

void Sampler::get_random_shortest_path()
{
    get_random_pair(s,t);
    find_shortest_path(s,t);
    if (layer.size() == 0) return;
    int32_t u,v,p;
    double sum;
    sum = 0;
    for (int32_t i = 0; i < layer.size(); ++i)
        sum += n_paths_from_s[layer[i]] * n_paths_from_t[layer[i]];
    sum *= rand_gen.get_real();
    for (int32_t i = 0; i < layer.size(); ++i)
    {
        sum -= n_paths_from_s[layer[i]] * n_paths_from_t[layer[i]];
        p = layer[i];
        if (sum < 0) break;
    }
    u = p;
    while (shortest_path_from_s[u] > 0)
    {
        path.push_back(u);
        layer.clear();
        sum = 0;
        neigh = g.get_adjlist_rev(u);
        for (int32_t i = 0; i < neigh.size(); ++i)
        {
            v = neigh[i];
            if (shortest_path_from_s[v] + 1 == shortest_path_from_s[u])
            {
                layer.push_back(v);
                sum += n_paths_from_s[v];
            }
        }
        sum *= rand_gen.get_real();
        for (int32_t i = 0; i < layer.size(); ++i)
        {
            sum -= n_paths_from_s[layer[i]];
            u = layer[i];
            if (sum < 0) break;
        }
        neigh.clear();
    }
//    path.push_back(s);
    reverse(path.begin(),path.end());
    u = p;
    while (shortest_path_from_t[u] > 1)
    {
        layer.clear();
        sum = 0;
        neigh = g.get_adjlist(u);
        for (int32_t i = 0; i < neigh.size(); ++i)
        {
            v = neigh[i];
            if (shortest_path_from_t[v] + 1 == shortest_path_from_t[u])
            {
                layer.push_back(v);
                sum += n_paths_from_t[v];
            }
        }
        sum *= rand_gen.get_real();
        for (int32_t i = 0; i < layer.size(); ++i)
        {
            sum -= n_paths_from_t[layer[i]];
            u = layer[i];
            if (sum < 0) break;
        }
        path.push_back(u);
        neigh.clear();
    }
//    if (p != t)
//        path.push_back(t);
    if (p == t) path.pop_back();
}

vector <int32_t> Sampler::bc_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
    n_fails = 0;
    while (1)
    {
        get_random_shortest_path();
        if (path.size() >= 1) break;
        ++n_fails;
    }
    zcover = 0;
    bool t1 = false, t2 = false;
    for (int i : path){
        if (seedMark[i]){
	    zcover = zcover | 2;
	    t1 = true;
	    if (t2)
		break;
	}
	if (cSeedMark[i]){
	    zcover = zcover | 1;
	    t2 = true;
	    if (t1)
		break;
	}
    }

    return path;
}

void Sampler::get_coverage_centrality_sample()
{
    vector <bool> in_set = vector <bool>(g.get_n_vertices()+1,false);
    get_random_pair(s,t);
    find_shortest_path(s,t);
    if (layer.size() == 0) return;
    in_set[s] = true;
    in_set[t] = true;
    queue <int32_t> q;
    int32_t u,v;
    //path.push_back(s);
    for (int32_t i = 0; i < layer.size(); ++i)
    {
        q.push(layer[i]);
        if (!in_set[layer[i]])
        {
            path.push_back(layer[i]);
            in_set[layer[i]] = true;
        }
    }
    while (!q.empty())
    {
        u = q.front();
        q.pop();
        neigh = g.get_adjlist_rev(u);
        for (int32_t i = 0; i < neigh.size(); ++i)
        {
            v = neigh[i];
            if ((shortest_path_from_s[v] != -1) && (shortest_path_from_s[v] + 1 == shortest_path_from_s[u]))
            {
                if (!in_set[v])
                {
                    path.push_back(v);
                    in_set[v] = true;
                    q.push(v);
                }
            }
        }
    }
    for (int32_t i = 0; i < layer.size(); ++i)
    {
        q.push(layer[i]);
        if (!in_set[layer[i]])
        {
            path.push_back(i);
            in_set[layer[i]] = true;
        }
    }
    while (!q.empty())
    {
        u = q.front();
        q.pop();
        neigh = g.get_adjlist(u);
        for (int32_t i = 0; i < neigh.size(); ++i)
        {
            v = neigh[i];
            if ((shortest_path_from_t[v] != -1) && (shortest_path_from_t[v] + 1 == shortest_path_from_t[u]))
            {
                if (!in_set[v])
                {
                    path.push_back(v);
                    in_set[v] = true;
                    q.push(v);
                }
            }
        }
    }
    //path.push_back(t);
}

vector <int32_t> Sampler::cc_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
    n_fails = 0;
    while (1)
    {
        get_coverage_centrality_sample();
        if (path.size() > 0) break;
        ++n_fails;
    }
    zcover = 0;
    bool t1 = false, t2 = false;
    for (int i : path){
        if (seedMark[i]){
            zcover = zcover | 2;
            t1 = true;
            if (t2)
                break;
        }
        if (cSeedMark[i]){
            zcover = zcover | 1;
            t2 = true;
            if (t1)
                break;
        }
    }
    return path;
}

bool Sampler::get_random_triangle()
{
    path.clear();
    int64_t r = rand_gen.get_max(g.get_probability(g.get_n_vertices()));
    int32_t left = 1, right = g.get_n_vertices(), mid;
    int32_t s,u,v;
    while (left <= right)
    {
        mid = (left + right) / 2;
        if (g.get_probability(mid) > r)
        {
            s = mid;
            right = mid - 1;
        }
        else left = mid + 1;
    }
    path.push_back(s);
    neigh = g.get_adjlist(s);
    u = rand_gen.get_max(neigh.size());
    do
    {
        v = rand_gen.get_max(neigh.size());
    } while (u == v);
    if (g.get_deg(u) > g.get_deg(v)) swap(u,v);
    u = neigh[u];v=neigh[v];   
	neigh.clear();
    neigh = g.get_adjlist(u);
    for (int32_t i = 0; i < neigh.size(); ++i)
        if (neigh[i] == v)
         {
             path.push_back(u);
             path.push_back(v);
             return true;
         }	    
    return false;
}

vector <int32_t> Sampler::triangle_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
    path.clear();
    n_fails = 0;
    while(!get_random_triangle())
    {
        ++n_fails;
    }
    zcover = 0;
    bool t1 = false, t2 = false;
    for (int i : path){
        if (seedMark[i]){
            zcover = zcover | 2;
            t1 = true;
            if (t2)
                break;
        }
        if (cSeedMark[i]){
            zcover = zcover | 1;
            t2 = true;
            if (t1)
                break;
        }
    }
    return path;
}

vector < pair<int32_t, int32_t> > Sampler::triangle_edge_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
    path.clear();
    n_fails = 0;
    while(!get_random_triangle())
    {
        ++n_fails;
    }
    vector < pair<int32_t, int32_t> > edges;
    edges.push_back(make_pair(path[0],path[1]));
    edges.push_back(make_pair(path[1],path[2]));
    edges.push_back(make_pair(path[2],path[0]));
    return edges;
}

void Sampler::init_k_path(int32_t init_k_path)
{
	k_path = init_k_path;
	explored = vector <bool>(g.get_n_vertices()+1,false);
}

vector <int32_t> Sampler::k_path_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
    n_fails = 0;
    path.clear();    
    int k_path_tmp = rand_gen.get_max(k_path) + 1;
    int32_t s;
    s = get_random_vertex();
    explored[s] = true;
    path.push_back(s);
    int32_t r;
    for (int32_t i = 0; i < k_path_tmp; ++i)
    {
        neigh = g.get_adjlist(s);
        r = neigh.size();
        for (int32_t j = 0; j < neigh.size(); ++j)
            if (explored[neigh[j]]) --r;
        if (r == 0) break;
        r = rand_gen.get_max(r) + 1;
        for (int32_t j = 0; j < neigh.size(); ++j)
        {
            if (explored[neigh[j]]) continue;
            --r;
            if (r == 0)
            {
                s = neigh[j];
                break;
            }
        }
        explored[s] = true;
        path.push_back(s);
        neigh.clear();
    }
	for (int i = 0; i < path.size(); ++i)
		explored[path[i]] = false;
    zcover = 0;
    bool t1 = false, t2 = false;
    for (int i : path){
        if (seedMark[i]){
            zcover = zcover | 2;
            t1 = true;
            if (t2)
                break;
        }
        if (cSeedMark[i]){
            zcover = zcover | 1;
            t2 = true;
            if (t1)
                break;
        }
    }
    return path;
}

vector <int32_t> Sampler::ds2_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
    path.clear();
    int s = rand_gen.get_max(g.get_n_vertices()) + 1;    
    neigh = g.get_adjlist(s);
    explored[s] = true;
    for (int v: neigh) 
    {
        path.push_back(v);
        explored[v]=true;
    }
    int l = path.size();
    int v;
    for (int32_t i = 0; i < l; ++i)
    {
        v = path[i];
        neigh = g.get_adjlist(v);
        for (int u: neigh)
            if (!explored[u])
            {
                path.push_back(u);
                explored[u] = false;
            }
    }
    path.push_back(s);
    for (int i = 0; i < path.size(); ++i)
        explored[path[i]] = false;
    n_fails = 0;
    zcover = 0;
    bool t1 = false, t2 = false;
    for (int i : path){
        if (seedMark[i]){
            zcover = zcover | 2;
            t1 = true;
            if (t2)
                break;
        }
        if (cSeedMark[i]){
            zcover = zcover | 1;
            t2 = true;
            if (t1)
                break;
        }
    }
    return path;
}
