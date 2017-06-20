#include "stdafx.h"
#include "Database.h"
#include "Ullman.h"
#include <cstring>
#include <string>
#include <iterator>
#include <strstream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include "timer.h"
#include <time.h> 
#include <ctime>
#include <fstream>
#include "GraphEditDistance.h"
#include <cmath>
#include <sys/timeb.h>		//Due to we have many applications(partitioning) in one second, we use microsecond as seed, if not portable we can use srand(time(NULL)), or do not use seed
///////////////Fisher-Yates shuffle/////////////////////////////////////
template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
//	srand(time(NULL));
	struct timeb tmb;
	ftime(&tmb);
	srand(tmb.millitm);		//Due to we have many applications(partitioning) in one second, we use microsecond as seed, if not portable we can use srand(time(NULL)), or do not use seed


	size_t left = std::distance(begin, end);
	while (num_random--) {
		bidiiter r = begin;
		std::advance(r, rand() % left);
		std::swap(*begin, *r);
		++begin;
		--left;
	}
	return begin;
}

template<class bidiiter>
bidiiter random_unique2(bidiiter begin, bidiiter end, size_t num_random) {
	srand(time(NULL));
	size_t left = std::distance(begin, end);
	while (num_random--) {
		bidiiter r = begin;
		std::advance(r, rand() % left);
		std::swap(*begin, *r);
		++begin;
		--left;
	}
	return begin;
}

namespace PSS{

	int Graph::buildEdge(){
		if ((*this).size() == 0) return -1;
		char buf[512];
		std::map <std::string, unsigned int> tmp;

		unsigned int id = 0;
		for (int from = 0; from < (int)size(); ++from) {
			for (Vertex::edge_iterator it = (*this)[from].edge.begin();
				it != (*this)[from].edge.end(); ++it)
			{
				//		cout << "label " << (*this)[from].label << " edge list size " << (*this)[from].edge.size() << endl;
				if (directed || from <= it->to)
					std::sprintf(buf, "%d %d %d", from, it->to, it->elabel);
				else
					std::sprintf(buf, "%d %d %d", it->to, from, it->elabel);

				// Assign unique id's for the edges.
				if (tmp.find(buf) == tmp.end()) {
					it->id = id;
					tmp[buf] = id;
					++id;
					//	cout << "id " <<id<< endl;
				}
				else {
					it->id = tmp[buf];
				}
			}
		}
		edge_size_ = id;
	}


template <class T, class Iterator>
void tokenize(const char *str, Iterator iterator)
{
	std::istrstream is(str, std::strlen(str));
	std::copy(std::istream_iterator <T>(is), std::istream_iterator <T>(), iterator);
}
void Graph::printGraph(std::string d){
	ofstream oss;
	oss.open(d, fstream::app | fstream::out);
	for (unsigned int i = 0; i < size(); i++){
		oss << "v" <<" "<< i <<" "<< (*this)[i].label << std::endl;
	}
	for (unsigned int i = 0; i < size(); i++){
		for (unsigned int j = 0; j < (*this)[i].edge.size(); j++){
			if ((*this)[i].edge[j].from<(*this)[i].edge[j].to)
			oss << "e " << (*this)[i].edge[j].from << " "  << (*this)[i].edge[j].to << " " << (*this)[i].edge[j].elabel << std::endl;
		}
	}
	oss << "t" << std::endl;
	oss.close();
}

void Graph::printGraph(std::string d, std::string code) {
	ofstream oss;
	oss.open(d, fstream::app | fstream::out);
	oss << code << endl;
	for (unsigned int i = 0; i < size(); i++){
		oss << "v" << " " << i << " " << (*this)[i].label << std::endl;
	}
	for (unsigned int i = 0; i < size(); i++){
		for (unsigned int j = 0; j < (*this)[i].edge.size(); j++){
			if ((*this)[i].edge[j].from<(*this)[i].edge[j].to)
				oss << "e " << (*this)[i].edge[j].from << " " << (*this)[i].edge[j].to << " " << (*this)[i].edge[j].elabel << std::endl;
		}
	}
	oss << "t" << std::endl;

	oss.close();
}

void Graph::buildProfile(){
	for (int j = 0; j < (*this).size(); j++){
		std::map <int, int>::iterator l;
		int c = (*this)[j].label;
		if (c == -2) continue;
		l = Label_statistic_v.find(c);
		if (l == Label_statistic_v.end()){
			std::pair<int, int> p;
			p.first = c;
			p.second = 1;
			Label_statistic_v.insert(p);
		}
		else{
			l->second++;
		}
		for (int k = 0; k < (*this)[j].edge.size(); k++){
			std::map <int, int>::iterator l;
			int c = (*this)[j].edge[k].elabel;
			if (c <= -2) continue;
			l = Label_statistic_e.find(c);
			if (l == Label_statistic_e.end()){
				std::pair<int, int> p;
				p.first = c;
				p.second = 1;
				Label_statistic_e.insert(p);
			}
			else{
				l->second++;
			}
		}
	}
}
void Graph::showGraph(){
	std::cout << "-------------------------" << std::endl;
	std::cout << "Graph:" << std::endl;
	for (unsigned int i = 0; i < size(); i++){
		for (unsigned int j = 0; j < (*this)[i].edge.size(); j++){
			std::cout << "e_id:" << (*this)[i].edge[j].id<< " " << (*this)[i].edge[j].from << "-" << (*this)[i].edge[j].to
				<< " uv_label:" << (*this)[i].label << "-" << (*this)[(*this)[i].edge[j].to].label << " e_lable:" << (*this)[i].edge[j].elabel << std::endl;
		}
	}
	std::cout << "-------------------------" << std::endl;
}


std::istream &Graph::read(std::istream &is)
{
	std::vector <std::string> result;
	char line[1024];

	clear();

	while (true) {

		unsigned int pos = is.tellg();  //std::cout << "pos " << pos << std::endl; //std::cin.get();
		if (!is.getline(line, 1024))
			break;
		result.clear();
		tokenize<std::string>(line, std::back_inserter(result));

		if (result.empty()) {
			// do nothing
		}
		else if (result[0] == "t") {
			if (!empty()) { // use as delimiter
//				std::cout << "t here" << std::endl;
//				is.seekg(pos, std::ios_base::beg);
				break;
			}
			else {
				/*
				* y = atoi (result[3].c_str());
				*/
			}
		}
		else if (result[0] == "v" && result.size() >= 3) {
			unsigned int id = atoi(result[1].c_str());//	std::cout << "v " << id << " - " << atoi(result[2].c_str()) << std::endl;
			this->resize(id + 1);
			(*this)[id].label = atoi(result[2].c_str());
		}
		else if (result[0] == "e" && result.size() >= 4) {
			int from = atoi(result[1].c_str());
			int to = atoi(result[2].c_str());
			int elabel = atoi(result[3].c_str());

			if ((int)size() <= from || (int)size() <= to) {
				std::cerr << "Format Error:  define vertex lists before edges" << std::endl;// cin.get();
				exit(-1);
			}

			(*this)[from].push(from, to, elabel); // cout << "from " << from << " to " << to << " elabel " << elabel << endl;
			if (directed == false)
			{
				(*this)[to].push(to, from, elabel);
				//		cout << "to " << to << " from " << from << " elabel " << elabel << endl;
			}
		}
	}
//	std::cout << "build edge here" << std::endl;
	buildEdge();

//	for (unsigned int i = 0; i < this->size(); i++){ Subgraph_Vertex_mapping.push_back(i); }
	return is;
}

bool in_vector(std::vector<int> V, int v){
	for (unsigned int i = 0; i < V.size(); i++) if (V[i] == v) return true;
	return false;
}

int vec_position(std::vector<int> V, int v){
	for (unsigned int i = 0; i < V.size(); i++) if (V[i] == v) return i;
}

std::vector < Graph > Database::Graph_partition_1(Graph &g, unsigned int k){
	if (k > g.size()) exit;
	std::vector<int> a(g.size());
	for (int i = 0; i<g.size(); ++i) a[i] = i;		
//	srand(time(NULL));
	random_unique(a.begin(), a.end(), k);
//	random_unique2(a.begin(), a.end(), a.size());
	/*
	for (int i = 0; i < a.size(); ++i) {
		std::cout << a[i] << ' ';
	}std::cout << '\n';
	*/
	std::vector<Graph> P;
	for (int i = 0; i<k; ++i) {
		Graph t_g; 
		t_g.push_back(g[a[i]]);
		t_g[0].edge.resize(0);
		P.push_back(t_g);
	}
	std::vector<std::vector<int>> P_node_mapping;			//Mapping of each partition's node to the node in the original graph
	P_node_mapping.resize(k);
	std::vector<int> b(g.size());
	for (int i = 0; i<g.size(); ++i) b[i] = -1;
	for (int i = 0; i < k; ++i) { 
		P_node_mapping[i].push_back(a[i]);
		b[a[i]] = i; 
	}


	int nodeCounter = 0; 
	int partionCounter=0;
	int t = 0;
	int loop = 0;
	while (nodeCounter < g.size() - k){
		int flag = 0;
		for (unsigned int i = 0; i < P_node_mapping[partionCounter].size(); i++){
			t = P_node_mapping[partionCounter][i];
			for (unsigned int j = 0; j < g[t].edge.size(); j++){
				if (b[(g[t].edge[j].to)] == -1) {
					Vertex v;
					v.label= g[(g[t].edge[j].to)].label;
					P[partionCounter].push_back(v);
					b[(g[t].edge[j].to)] = partionCounter;
					P_node_mapping[partionCounter].push_back(g[t].edge[j].to);
					nodeCounter++;
					flag = 1;
					break;
				}
			}
			if (flag == 1) break;
		}

		partionCounter++; 
		if (flag == 1) loop = 0;
		if ((flag == 0) && (partionCounter > P_node_mapping.size()-1)) { loop++; }
		if ((flag == 0) && (loop == 2)) {
			//	Some data graphs are disconnected, this makes the graph canonicalization tool problematic.
			//	Here when we meet a disconnected graph, we will first find its largest partitions,
			//	then we report the graph induced by nodes in each partition as a partition.
			//	As a result, some nodes and edges are discarded, however the pigeon hold principle still holds the query processing is still correct.
			for (unsigned int i = 0; i < P_node_mapping.size(); i++){
				for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
					for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
						if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
							Edge e = g[P_node_mapping[i][j]].edge[k];
							e.from = vec_position(P_node_mapping[i], e.from);
							if (e.from != j) std::cout << "error" << std::endl;
							e.to = vec_position(P_node_mapping[i], e.to);
							if (e.from == e.to) std::cout << "error2" << std::endl;
							P[i][j].edge.push_back(e);
						}
					}
				}
			} 
			for (unsigned int i = 0; i < P.size(); i++){
				P[i].buildEdge();
			}
			return P;
		}
		partionCounter = partionCounter % k; 
//		std::cout << nodeCounter << " " << flag << "|" << loop << "  " << partionCounter << "|" << P_node_mapping.size() << '\r';
	} 
	/*
	for (int i = 0; i < b.size(); ++i) {
		std::cout << b[i] << ' ';
	}std::cout << '\n';
	
	for (int i = 0; i<g[0].edge.size(); ++i) {
		std::cout << g[0].edge[i].from << ' ' << g[0].edge[i].to << std::endl;
	}std::cout << '\n';
	std::cout << "here4" << std::endl; std::cin.get();
	*/
	std::map<std::pair<int,int>, int> halfEdge_map;
		
		for (unsigned int i = 0; i < P_node_mapping.size(); i++){
			for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
//				std::cout << P_node_mapping[i][j] << '-' << g[P_node_mapping[i][j]].edge.size()<<"  ";
//				P[i][j].edge.resize(0);
				for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
					if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){ 
						Edge e = g[P_node_mapping[i][j]].edge[k];
						e.from = vec_position(P_node_mapping[i],e.from);
						if (e.from != j) std::cout << "error" << std::endl;
						e.to = vec_position(P_node_mapping[i], e.to);
						if (e.from == e.to) std::cout << "error2" << std::endl;
						P[i][j].edge.push_back(e); 
					}
					else{
////////////Here we assign half edges to each partition////////////////////////////////////////////
						int t = -1;
						for (unsigned int l = 0; l < P_node_mapping.size(); l++){
							if (in_vector(P_node_mapping[l], g[P_node_mapping[i][j]].edge[k].to)){ if (i == l) std::cout << "errpr3" << std::endl; t = l; }
						}
						std::pair<int,int> a, b;
							a.first = g[P_node_mapping[i][j]].edge[k].from;
							a.second = g[P_node_mapping[i][j]].edge[k].to;
							b.first = g[P_node_mapping[i][j]].edge[k].to;
							b.second = g[P_node_mapping[i][j]].edge[k].from;
					
							if (rand() % 2 == 0){ 							
							halfEdge_map.insert(std::pair<std::pair<int,int>, int>(a, i));
							halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, i));
						}
						else{
							halfEdge_map.insert(std::pair<std::pair<int, int>, int>(a, t));
							halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, t));
						}
					}
				}
			}//std::cout << '\n';
		}
		/*
	
		std::cout << '\n';
		for (unsigned int i = 0; i < P_node_mapping.size(); i++){
			for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
				std::cout << P[i][j].edge.size() << ' ';
				
			}std::cout << '\n';
		}
		std::cout << '\n';
*/
		/*
		std::cout << '\n';
		for (unsigned int i = 0; i < P_node_mapping.size(); i++){
			for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
//				P[i][j].showEdge();
				std::cout << P_node_mapping[i][j] << "|" << g[P_node_mapping[i][j]].label << "|" << P[i][j].label << ' ';
			}std::cout << '\n';
		}
	*/
		for (unsigned int i = 0; i < P_node_mapping.size(); i++){
			Vertex v;
			v.label = -2;		//-2 means it is a dummy node all half edge points to it and it is always the last vertex in a graph
			v.edge.resize(0);
			P[i].push_back(v);
			for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
				int elabel = -2;
				for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
					if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
					}
					else{
						////////////Here we assign each half edges to each vertex////////////////////////////////////////////
						std::pair<int, int> p_t;
						p_t.first = g[P_node_mapping[i][j]].edge[k].from;
						p_t.second = g[P_node_mapping[i][j]].edge[k].to;
						if (halfEdge_map.find(p_t)->second == i){
							Edge e1,e2;
							e1.from = j;
							e1.to = P[i].size()-1;
							e1.elabel = elabel;
							e2.to = j;
							e2.from = P[i].size() - 1;
							e2.elabel = elabel;
							elabel--;
							/////////Here we add the half edge to the dummy node
							P[i][j].edge.push_back(e1);
							P[i][P[i].size() - 1].edge.push_back(e2);
						}
					}
				}
			}
		}
		for (unsigned int i = 0; i < P.size(); i++){
			if (P[i][P[i].size() - 1].edge.size() == 0)	P[i].pop_back();		// If a partition has no half edge at all we will delete the dummy node
		}
		/*
		for (unsigned int i = 0; i < P.size(); i++){
			for (unsigned int j = 0; j < P[i].size(); j++){
				std::cout << P[i][j].edge.size() << ' ';

			}std::cout << '\n';
		}
		std::cout << '\n';
*/
		for (unsigned int i = 0; i < P.size(); i++){
			P[i].buildEdge();
		}
		/*
		for (int i = 0; i < b.size(); ++i) {
			std::cout << b[i] << ' ';
		}std::cout << '\n';
*/
		vector<int> zzzz; int f = 1;
		for (unsigned int i = 0; i < P_node_mapping.size(); i++){
			for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
				//		cout << P_node_mapping[i][j] << " ";// "-" << g[P_node_mapping[i][j]].label

				for (int n = 0; n < zzzz.size(); n++){
					if (zzzz[n] == P_node_mapping[i][j]) { std::cout << "Wrong!!!!!!!!!" << endl; cin.get(); f = 0; }

				}zzzz.push_back(P_node_mapping[i][j]);
			}
		}
		return P;
}

std::vector < Graph > Database::Graph_partition_2(Graph &g, unsigned int k){
	if (k > g.size()) exit;
	std::vector<int> a(g.size());
	for (int i = 0; i<g.size(); ++i) a[i] = i;
	random_unique(a.begin(), a.end(), k);
	std::vector<Graph> P;
	for (int i = 0; i<k; ++i) {
		Graph t_g;
		t_g.push_back(g[a[i]]);
		t_g[0].edge.resize(0);
		P.push_back(t_g);
	}
	std::vector<std::vector<int>> P_node_mapping;			//Mapping of each partition's node to the node in the original graph
	P_node_mapping.resize(k);
	std::vector<int> b(g.size());
	for (int i = 0; i<g.size(); ++i) b[i] = -1;
	for (int i = 0; i < k; ++i) {
		P_node_mapping[i].push_back(a[i]);
		b[a[i]] = i;
	}


	int nodeCounter = 0;
	int partionCounter = 0;
	int t = 0;
	int loop = 0;
	while (nodeCounter < g.size() - k){
		int flag = 0;
		for (unsigned int i = 0; i < P_node_mapping[partionCounter].size(); i++){
			t = P_node_mapping[partionCounter][P_node_mapping[partionCounter].size()-i-1];
			for (unsigned int j = 0; j < g[t].edge.size(); j++){
				if (b[(g[t].edge[j].to)] == -1) {
					Vertex v;
					v.label = g[(g[t].edge[j].to)].label;
					P[partionCounter].push_back(v);
					b[(g[t].edge[j].to)] = partionCounter;
					P_node_mapping[partionCounter].push_back(g[t].edge[j].to);
					nodeCounter++;
					flag = 1;
					break;
				}
			}
			if (flag == 1) break;
		}

		partionCounter++;
		if (flag == 1) loop = 0;
		if ((flag == 0) && (partionCounter > P_node_mapping.size() - 1)) { loop++; }
		if ((flag == 0) && (loop == 2)) {
			//	Some data graphs are disconnected, this makes the graph canonicalization tool problematic.
			//	Here when we meet a disconnected graph, we will first find its largest partitions,
			//	then we report the graph induced by nodes in each partition as a partition.
			//	As a result, some nodes and edges are discarded, however the pigeon hold principle still holds the query processing is still correct.
			for (unsigned int i = 0; i < P_node_mapping.size(); i++){
				for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
					for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
						if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
							Edge e = g[P_node_mapping[i][j]].edge[k];
							e.from = vec_position(P_node_mapping[i], e.from);
							if (e.from != j) std::cout << "error" << std::endl;
							e.to = vec_position(P_node_mapping[i], e.to);
							if (e.from == e.to) std::cout << "error2" << std::endl;
							P[i][j].edge.push_back(e);
						}
					}
				}
			}
			for (unsigned int i = 0; i < P.size(); i++){
				P[i].buildEdge();
			}
			return P;
		}
		partionCounter = partionCounter % k;
		//		std::cout << nodeCounter << " " << flag << "|" << loop << "  " << partionCounter << "|" << P_node_mapping.size() << '\r';
	}

	std::map<std::pair<int, int>, int> halfEdge_map;

	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
			//				std::cout << P_node_mapping[i][j] << '-' << g[P_node_mapping[i][j]].edge.size()<<"  ";
			//				P[i][j].edge.resize(0);
			for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
				if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
					Edge e = g[P_node_mapping[i][j]].edge[k];
					e.from = vec_position(P_node_mapping[i], e.from);
					if (e.from != j) std::cout << "error" << std::endl;
					e.to = vec_position(P_node_mapping[i], e.to);
					if (e.from == e.to) std::cout << "error2" << std::endl;
					P[i][j].edge.push_back(e);
				}
				else{
					////////////Here we assign half edges to each partition////////////////////////////////////////////
					int t = -1;
					for (unsigned int l = 0; l < P_node_mapping.size(); l++){
						if (in_vector(P_node_mapping[l], g[P_node_mapping[i][j]].edge[k].to)){ if (i == l) std::cout << "errpr3" << std::endl; t = l; }
					}
					std::pair<int, int> a, b;
					a.first = g[P_node_mapping[i][j]].edge[k].from;
					a.second = g[P_node_mapping[i][j]].edge[k].to;
					b.first = g[P_node_mapping[i][j]].edge[k].to;
					b.second = g[P_node_mapping[i][j]].edge[k].from;

					if (rand() % 2 == 0){
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(a, i));
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, i));
					}
					else{
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(a, t));
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, t));
					}
				}
			}
		}
	}

	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		Vertex v;
		v.label = -2;		//-2 means it is a dummy node all half edge points to it and it is always the last vertex in a graph
		v.edge.resize(0);
		P[i].push_back(v);
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
			int elabel=-2;
			for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
				if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
				}
				else{
					////////////Here we assign each half edges to each vertex////////////////////////////////////////////
					std::pair<int, int> p_t;
					p_t.first = g[P_node_mapping[i][j]].edge[k].from;
					p_t.second = g[P_node_mapping[i][j]].edge[k].to;
					if (halfEdge_map.find(p_t)->second == i){
						Edge e1, e2;
						e1.from = j;
						e1.to = P[i].size() - 1;
						e1.elabel = elabel;
						e2.to = j;
						e2.from = P[i].size() - 1;
						e2.elabel = elabel;
						elabel--;
						/////////Here we add the half edge to the dummy node
						P[i][j].edge.push_back(e1);
						P[i][P[i].size() - 1].edge.push_back(e2);
					}
				}
			}
		}
	}
	for (unsigned int i = 0; i < P.size(); i++){
		if (P[i][P[i].size() - 1].edge.size() == 0)	P[i].pop_back();		// If a partition has no half edge at all we will delete the dummy node
	}
	for (unsigned int i = 0; i < P.size(); i++){
		P[i].buildEdge();
	}
	return P;
}

std::vector < Graph > Database::Graph_partition_3(Graph &g, unsigned int k){
	if (k > g.size()) exit;
	std::vector<int> a(g.size());
	for (int i = 0; i<g.size(); ++i) a[i] = i;
	random_unique(a.begin(), a.end(), k);					// Randomly assign 
	
	for (int i = 0; i < k; i++){							// Choose a least frequent label node from the neighbors of a[i] as well as a[i] itself
		int temp_node_min_fre = Label_statistic_v.find(g[a[i]].label)->second;
		int min = a[i];
		for (int j = 0; j < g[a[i]].edge.size(); j++){
			if (Label_statistic_v.find(g[(g[a[i]].edge[j].to)].label)->second < temp_node_min_fre){
				temp_node_min_fre = Label_statistic_v.find(g[(g[a[i]].edge[j].to)].label)->second;
				int flag = 1;
				for (int b = 0; b < k; b++){
					if (a[b] == g[a[i]].edge[j].to) flag = 0;
				}
				if (flag) min = g[a[i]].edge[j].to;
			}
		}
		a[i] = min;
	}

	std::vector<Graph> P;									// P is a vector of k graphs
	for (int i = 0; i<k; ++i) {
		Graph t_g;
		t_g.push_back(g[a[i]]);
		t_g[0].edge.resize(0);
		P.push_back(t_g);
	}
	std::vector<std::vector<int>> P_node_mapping;			// Mapping of each partition's node to the node in the original graph
	P_node_mapping.resize(k);
	std::vector<int> b(g.size());							// Record the use of each node in G
	for (int i = 0; i<g.size(); ++i) b[i] = -1;
	for (int i = 0; i < k; ++i) {
		P_node_mapping[i].push_back(a[i]);
		b[a[i]] = i;
	}

	std::vector<double> SelectivityValue;					// 
	SelectivityValue.resize(k);
	for (int i = 0; i < k; i++){	
		SelectivityValue[i] = Label_statistic_v.find(g[a[i]].label)->second;
	}
	std::vector<double> SumFreE(k,0);					// Sum of the frequence of all the edges in current partition		
	std::vector<double> SumFreV;								// Sum of the frequence of all the nodes in current partition
	SumFreV = SelectivityValue;
	std::vector<int> VSize(k,1);						// Size is #nodes in current partition	
	std::vector<int> ESize(k, 0);						// Size is #edge in current partition	
	std::vector<double> NewSelectivityValue=SelectivityValue;			// Record each partition's optimal new SelectivityValue

	int nodeCounter = 0;
	
	int t = 0;
	
	while (nodeCounter < g.size() - k){
		std::vector<double> OptimalGrowth(k, DBL_MAX);					// OptimalGrowth is a value for each partition when deciding which partition and which node to growth that has best result
		std::vector<int> OptimalGrowthNode;						// 
		OptimalGrowthNode.resize(k);
		int flag = 0;
		for (int partionCounter = 0; partionCounter < k;partionCounter++){		
			for (unsigned int i = 0; i < P_node_mapping[partionCounter].size(); i++){
				t = P_node_mapping[partionCounter][i];
				for (unsigned int j = 0; j < g[t].edge.size(); j++){
					if (b[(g[t].edge[j].to)] == -1) {
						flag = 1;
						int Vfre=Label_statistic_v.find(g[(g[t].edge[j].to)].label)->second;
						int Efre = 0;
						int edgeCounter = 0;

						for (int k = 0; k < g[(g[t].edge[j].to)].edge.size(); k++){
							if (b[g[(g[t].edge[j].to)].edge[k].to] == partionCounter){
								Efre += Label_statistic_e.find(g[(g[t].edge[j].to)].edge[k].elabel)->second;
								edgeCounter++;
							}					
						}

						double Select = 0,temp;
						temp = ((SumFreV[partionCounter] + Vfre) / (VSize[partionCounter] + 1) + (SumFreE[partionCounter] + Efre) / (ESize[partionCounter] + edgeCounter)) / pow((VSize[partionCounter] + 1 + ESize[partionCounter] + edgeCounter),2);
						 Select = temp-SelectivityValue[partionCounter];
					
						if (Select < OptimalGrowth[partionCounter])
						{
				//			std::cout << "Better Sel:" << Select << '\n';
				//			std::cout << "Better Node:" << g[t].edge[j].to << '\n';
							OptimalGrowth[partionCounter] = Select;
							OptimalGrowthNode[partionCounter] = g[t].edge[j].to;
							NewSelectivityValue[partionCounter] = temp;
						}
					}
				}			
			}
		}
		if (flag == 0){
			//	Some data graphs are disconnected, this makes the graph canonicalization tool problematic.
			//	Here when we meet a disconnected graph, we will first find its largest partitions,
			//	then we report the graph induced by nodes in each partition as a partition.
			//	As a result, some nodes and edges are discarded, however the pigeon hold principle still holds the query processing is still correct.
			//			for (unsigned int i = 0; i < P.size(); i++){
			//				P[i].buildEdge();
			//			}
			break;
		}
		int bestNode = 0, bestPartition=0;
		double bestValue = DBL_MAX;
		for (int i = 0; i<k; i++){
			if (OptimalGrowth[i] < bestValue) {
				bestValue = OptimalGrowth[i];
				bestNode = OptimalGrowthNode[i];
				bestPartition = i;
//				cout << "OptimalGrowth[i]:" << OptimalGrowth[i] << endl;
//				cout << "OptimalGrowthNode[i]:" << OptimalGrowthNode[i] << endl;
			}
		}
//		cout << "Best p:" << bestPartition << endl;
		Vertex v;
		v.label = g[bestNode].label;
		int edgeCounter = 0;
		int Efre = 0;
		for (int i = 0; i<g[bestNode].edge.size(); i++){
			if (b[g[bestNode].edge[i].to] == bestPartition){
//				Edge e;
//				e = g[bestNode].edge[i];
//				e.from = vec_position(P_node_mapping[bestPartition], e.from);
//				e.to = vec_position(P_node_mapping[bestPartition], e.to);
//				v.edge.push_back(e);
				edgeCounter++;
				Efre += Label_statistic_e.find(g[bestNode].edge[i].elabel)->second;
			}
		}
		VSize[bestPartition] += 1;				// 1 new node
		ESize[bestPartition] += edgeCounter ;  // #new edges 
		int Vfre = Label_statistic_v.find(g[bestNode].label)->second;
		P[bestPartition].push_back(v);
		SumFreV[bestPartition] += Vfre;
		SumFreE[bestPartition] += Efre;
		P_node_mapping[bestPartition].push_back(bestNode);
		b[bestNode] = bestPartition;
		nodeCounter++;
		SelectivityValue[bestPartition] = NewSelectivityValue[bestPartition];
	}
	/*
	for (int i = 0; i < b.size(); ++i) {
	std::cout << b[i] << ' ';
	}std::cout << '\n';
*/

/*
	for (int i = 0; i<g[0].edge.size(); ++i) {
	std::cout << g[0].edge[i].from << ' ' << g[0].edge[i].to << std::endl;
	}std::cout << '\n';
*/	
	std::map<std::pair<int, int>, int> halfEdge_map;
	
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
			//				std::cout << P_node_mapping[i][j] << '-' << g[P_node_mapping[i][j]].edge.size()<<"  ";
			//				
			P[i][j].edge.resize(0);
			for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
				if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
					Edge e = g[P_node_mapping[i][j]].edge[k];
					e.from = vec_position(P_node_mapping[i], e.from);
					if (e.from != j) std::cout << "error" << std::endl;
					e.to = vec_position(P_node_mapping[i], e.to);
					if (e.from == e.to) std::cout << "error2" << std::endl;
					P[i][j].edge.push_back(e);
				}
				else{
					////////////Here we assign half edges to each partition////////////////////////////////////////////
					int t = -1;
					for (unsigned int l = 0; l < P_node_mapping.size(); l++){
						if (in_vector(P_node_mapping[l], g[P_node_mapping[i][j]].edge[k].to)){ if (i == l) std::cout << "errpr3" << std::endl; t = l; }
					}
					std::pair<int, int> a, b;
					a.first = g[P_node_mapping[i][j]].edge[k].from;
					a.second = g[P_node_mapping[i][j]].edge[k].to;
					b.first = g[P_node_mapping[i][j]].edge[k].to;
					b.second = g[P_node_mapping[i][j]].edge[k].from;

					if (rand() % 2 == 0){
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(a, i));
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, i));
					}
					else{
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(a, t));
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, t));
					}
				}
			}
		}//std::cout << '\n';
	}
	/*

	std::cout << '\n';
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
	for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
	std::cout << P[i][j].edge.size() << ' ';

	}std::cout << '\n';
	}
	std::cout << '\n';
	*/
	/*
	std::cout << '\n';
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
	for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
	//				P[i][j].showEdge();
	std::cout << P_node_mapping[i][j] << "|" << g[P_node_mapping[i][j]].label << "|" << P[i][j].label << ' ';
	}std::cout << '\n';
	}
	*/
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		Vertex v;
		v.label = -2;		//-2 means it is a dummy node all half edge points to it and it is always the last vertex in a graph
		v.edge.resize(0);
		P[i].push_back(v);
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
			int elabel=-2;
			for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
				if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
				}
				else{
					////////////Here we assign each half edges to each vertex////////////////////////////////////////////
					std::pair<int, int> p_t;
					p_t.first = g[P_node_mapping[i][j]].edge[k].from;
					p_t.second = g[P_node_mapping[i][j]].edge[k].to;
					if (halfEdge_map.find(p_t)->second == i){
						Edge e1, e2;
						e1.from = j;
						e1.to = P[i].size() - 1;
						e1.elabel = elabel;
						e2.to = j;
						e2.from = P[i].size() - 1;
						e2.elabel = elabel;
						elabel--;
						/////////Here we add the half edge to the dummy node
						P[i][j].edge.push_back(e1);
						P[i][P[i].size() - 1].edge.push_back(e2);
					}
				}
			}
		}
	}
	for (unsigned int i = 0; i < P.size(); i++){
		if (P[i][P[i].size() - 1].edge.size() == 0)	P[i].pop_back();		// If a partition has no half edge at all we will delete the dummy node
	}
	/*
	for (unsigned int i = 0; i < P.size(); i++){
	for (unsigned int j = 0; j < P[i].size(); j++){
	std::cout << P[i][j].edge.size() << ' ';

	}std::cout << '\n';
	}
	std::cout << '\n';
	*/
	for (unsigned int i = 0; i < P.size(); i++){
		P[i].buildEdge();
	}
	/*
	for (int i = 0; i < b.size(); ++i) {
	std::cout << b[i] << ' ';
	}std::cout << '\n';
	

	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		cout << "P" << i << "  Sel:"<<SelectivityValue[i]<<"   ";
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
			cout << P_node_mapping[i][j] <<  " ";// "-" << g[P_node_mapping[i][j]].label
		}
		cout << endl;
	}
	cin.get();
*/
	vector<int> zzzz; int f = 1;
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
	//		cout << P_node_mapping[i][j] << " ";// "-" << g[P_node_mapping[i][j]].label

			for (int n = 0; n < zzzz.size(); n++){
				if (zzzz[n] == P_node_mapping[i][j]) { std::cout << "Wrong!!!!!!!!!" << endl; cin.get(); f = 0; }
				
			}zzzz.push_back(P_node_mapping[i][j]);
		}
	}

	return P;
}
std::vector < Graph > Database::Graph_partition_4(Graph &g, unsigned int k){
//	cout << endl << "Graph ID:" << g.ID << "\r";
	int maxid = 0;
	for (int i = 0; i<g.size(); i++){
		for (int j = 0; j<g[i].edge.size(); j++){
	//		g[i].edge[j].show();
			if (g[i].edge[j].id>maxid) maxid = g[i].edge[j].id;			//Finding the number of edges in a graph
		}
	}
//	cout << "max " << maxid << endl;
//	cin.get();
	if (k > g.size()) { std::cout << "Graph size than number of partitions" <<g.ID<< endl; exit; }
	std::vector<int> a(g.size());
	for (int i = 0; i<g.size(); ++i) a[i] = i;
	random_unique(a.begin(), a.end(), k);					// Randomly assign 

	for (int i = 0; i < k; i++){							// Choose a least frequent label node from the neighbors of a[i] as well as a[i] itself
		int temp_node_min_fre = Label_statistic_v.find(g[a[i]].label)->second;
		int min = a[i];
		for (int j = 0; j < g[a[i]].edge.size(); j++){
			if (Label_statistic_v.find(g[(g[a[i]].edge[j].to)].label)->second < temp_node_min_fre){
				temp_node_min_fre = Label_statistic_v.find(g[(g[a[i]].edge[j].to)].label)->second;
				int flag = 1;
				for (int b = 0; b < k; b++){
					if (a[b] == g[a[i]].edge[j].to) flag = 0;
				}
				if (flag) min = g[a[i]].edge[j].to;
			}
		}
		a[i] = min;
	}
	/*
	for (int i = 0; i < k; i++){
		std::cout << a[i] << " ";
	}
	std::cout << std::endl;
	cin.get();
*/

	std::vector<Graph> P;									// P is a vector of k graphs
	for (int i = 0; i<k; ++i) {
		Graph t_g;
		t_g.push_back(g[a[i]]);
		t_g[0].edge.resize(0);
		P.push_back(t_g);
	}
	std::vector<std::vector<int>> P_node_mapping;			// Mapping of each partition's node to the node in the original graph
	P_node_mapping.resize(k);
	std::vector<int> b(g.size());							// Record the use of each node in G
	for (int i = 0; i<g.size(); ++i) b[i] = -1;
	for (int i = 0; i < k; ++i) {
		P_node_mapping[i].push_back(a[i]);
		b[a[i]] = i;
	}
	/*
	for (int i = 0; i < g.size(); i++){
	std::cout << b[i] << " ";
	}
	std::cout << std::endl;
	cin.get();
*/


	std::vector<double> SelectivityValue;					// 
	SelectivityValue.resize(k);
	for (int i = 0; i < k; i++){
		SelectivityValue[i] = 1.0/Label_statistic_v.find(g[a[i]].label)->second;
	}

	std::vector<double> SumFreE(k, 0);					// Sum of the frequence of all the edges in current partition		
	std::vector<double> SumFreV;								// Sum of the frequence of all the nodes in current partition
	SumFreV = SelectivityValue;

	std::vector<int> VSize(k, 1);						// Size is #nodes in current partition	
	std::vector<int> ESize(k, 0);						// Size is #edge in current partition	
	std::vector<double> NewSelectivityValue = SelectivityValue;			// Record each partition's optimal new SelectivityValue

	int nodeCounter = 0;

	while (nodeCounter < g.size() - k){
		std::vector<double> OptimalGrowth(k, -99);					// OptimalGrowth is a value for each partition when deciding which partition and which node to growth that has best result
		std::vector<int> OptimalGrowthNode(k, -1);						// 
		/*
		for (int i = 0; i < OptimalGrowthNode.size(); i++){
			std::cout << OptimalGrowthNode[i] << " ";
		}
		std::cout << endl;
		cin.get();
*/
		int flag = 0;
		for (int partionCounter = 0; partionCounter < k; partionCounter++){
			for (unsigned int i = 0; i < P_node_mapping[partionCounter].size(); i++){
				int t = P_node_mapping[partionCounter][i];
				for (unsigned int j = 0; j < g[t].edge.size(); j++){
					if (b[(g[t].edge[j].to)] == -1) {
						flag = 1;
						double Vfre = 1.0/Label_statistic_v.find(g[(g[t].edge[j].to)].label)->second;						
						double Efre = 0;
						int edgeCounter = 0;
						for (int k = 0; k < g[(g[t].edge[j].to)].edge.size(); k++){
							if (b[g[(g[t].edge[j].to)].edge[k].to] == partionCounter){
								Efre += 1.0/Label_statistic_e.find(g[(g[t].edge[j].to)].edge[k].elabel)->second;
								edgeCounter++;
							}
						}
						
						double Select = 0, temp;
						temp = ((SumFreV[partionCounter] + Vfre) / (VSize[partionCounter] + 1) + (SumFreE[partionCounter] + Efre) / (ESize[partionCounter] + edgeCounter)) * pow((double)(VSize[partionCounter] + 1 + ESize[partionCounter] + edgeCounter) / (double)(maxid + 1 + g.size()), 0.5);
			//			temp = ((SumFreV[partionCounter] + Vfre) + (SumFreE[partionCounter] + Efre)) * pow((double)(VSize[partionCounter] + 1 + ESize[partionCounter] + edgeCounter) / (double)(maxid + 1 + g.size()), 0.5);
						if (temp == 0)
						{
							cout << "Error temp" << endl; 
							cout << ((SumFreV[partionCounter] + Vfre) + (SumFreE[partionCounter] + Efre)) << " * " << (double)(VSize[partionCounter] + 1 + ESize[partionCounter] + edgeCounter) / (double)(maxid + 1 + g.size()) << endl;
							cout << (VSize[partionCounter] + 1 + ESize[partionCounter] + edgeCounter) <<" / "<< (maxid + 1 + g.size()) << endl;
							cin.get();
						}
						Select = temp - SelectivityValue[partionCounter];
		//				if (Select>0)
		//				cout << partionCounter << "     Select:" << Select << " T:" << temp << "   Optimal: " << OptimalGrowth[partionCounter] << "SelectivityValue[partionCounter]" << SelectivityValue[partionCounter] << endl;
						if (Select > OptimalGrowth[partionCounter])
						{
							//			std::cout << "Better Sel:" << Select << '\n';
							//			std::cout << "Better Node:" << g[t].edge[j].to << '\n';
							OptimalGrowth[partionCounter] = Select;
							OptimalGrowthNode[partionCounter] = g[t].edge[j].to;
							NewSelectivityValue[partionCounter] = temp;
						}
					}
				}
			}
		}
		if (flag == 0){
			//	Some data graphs are disconnected, this makes the graph canonicalization tool problematic.
			//	Here when we meet a disconnected graph, we will first find its largest partitions,
			//	then we report the graph induced by nodes in each partition as a partition.
			//	As a result, some nodes and edges are discarded, however the pigeon hold principle still holds the query processing is still correct.
			//			for (unsigned int i = 0; i < P.size(); i++){
			//				P[i].buildEdge();
			//			}
			break;
		}


		int bestNode = -1, bestPartition = -1;
		double bestValue = -99;
		
		for (int i = 0; i<k; i++){
			
			if (OptimalGrowth[i] > bestValue) {
				
				bestValue = OptimalGrowth[i];
				bestNode = OptimalGrowthNode[i];
				bestPartition = i;
				//				cout << "OptimalGrowth[i]:" << OptimalGrowth[i] << endl;
				//				cout << "OptimalGrowthNode[i]:" << OptimalGrowthNode[i] << endl;
			}
		}
		//		cout << "Best p:" << bestPartition << endl;
		Vertex v;
		v.label = g[bestNode].label;
		int edgeCounter = 0;
		double Efre = 0;
		for (int i = 0; i<g[bestNode].edge.size(); i++){
			if (b[g[bestNode].edge[i].to] == bestPartition){
				//				Edge e;
				//				e = g[bestNode].edge[i];
				//				e.from = vec_position(P_node_mapping[bestPartition], e.from);
				//				e.to = vec_position(P_node_mapping[bestPartition], e.to);
				//				v.edge.push_back(e);
				edgeCounter++;
				Efre += 1.0/Label_statistic_e.find(g[bestNode].edge[i].elabel)->second;
			}
		}
		VSize[bestPartition] += 1;				// 1 new node
		ESize[bestPartition] += edgeCounter;	// #new edges 
		double Vfre = 1.0/Label_statistic_v.find(g[bestNode].label)->second;
		P[bestPartition].push_back(v);
		SumFreV[bestPartition] += Vfre;
		SumFreE[bestPartition] += Efre;
//		cout << SumFreV[bestPartition] << "    " << SumFreE[bestPartition] << endl;  cin.get();
		P_node_mapping[bestPartition].push_back(bestNode);
//		cout << bestPartition << "    " << bestNode << endl;  
		b[bestNode] = bestPartition;
		nodeCounter++;
		SelectivityValue[bestPartition] = NewSelectivityValue[bestPartition];
	}
	/*
	for (int i = 0; i < b.size(); ++i) {
	std::cout << b[i] << ' ';
	}std::cout << '\n';
	*/

	/*
	for (int i = 0; i<g[0].edge.size(); ++i) {
	std::cout << g[0].edge[i].from << ' ' << g[0].edge[i].to << std::endl;
	}std::cout << '\n';
	*/
	std::map<std::pair<int, int>, int> halfEdge_map;

	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
			//				std::cout << P_node_mapping[i][j] << '-' << g[P_node_mapping[i][j]].edge.size()<<"  ";
			//				
			P[i][j].edge.resize(0);
			for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
				if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
					Edge e = g[P_node_mapping[i][j]].edge[k];
					e.from = vec_position(P_node_mapping[i], e.from);
					if (e.from != j) std::cout << "error" << std::endl;
					e.to = vec_position(P_node_mapping[i], e.to);
					if (e.from == e.to) std::cout << "error2" << std::endl;
					P[i][j].edge.push_back(e);
				}
				else{
					////////////Here we assign half edges to each partition////////////////////////////////////////////
					int t = -1;
					for (unsigned int l = 0; l < P_node_mapping.size(); l++){
						if (in_vector(P_node_mapping[l], g[P_node_mapping[i][j]].edge[k].to)){ if (i == l) std::cout << "error3" << std::endl; t = l; }
					}
					std::pair<int, int> a, b;
					a.first = g[P_node_mapping[i][j]].edge[k].from;
					a.second = g[P_node_mapping[i][j]].edge[k].to;
					b.first = g[P_node_mapping[i][j]].edge[k].to;
					b.second = g[P_node_mapping[i][j]].edge[k].from;

					if (rand() % 2 == 0){
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(a, i));
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, i));
					}
					else{
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(a, t));
						halfEdge_map.insert(std::pair<std::pair<int, int>, int>(b, t));
					}
				}
			}
		}//std::cout << '\n';
	}
	/*

	std::cout << '\n';
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
	for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
	std::cout << P[i][j].edge.size() << ' ';

	}std::cout << '\n';
	}
	std::cout << '\n';
	*/
	/*
	std::cout << '\n';
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
	for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
	//				P[i][j].showEdge();
	std::cout << P_node_mapping[i][j] << "|" << g[P_node_mapping[i][j]].label << "|" << P[i][j].label << ' ';
	}std::cout << '\n';
	}
	*/
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		Vertex v;
		v.label = -2;		//-2 means it is a dummy node all half edge points to it and it is always the last vertex in a graph
		v.edge.resize(0);
		P[i].push_back(v);
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
			int elabel = -2;
			for (unsigned int k = 0; k < g[P_node_mapping[i][j]].edge.size(); k++){
				if (in_vector(P_node_mapping[i], g[P_node_mapping[i][j]].edge[k].to)){
				}
				else{
					////////////Here we assign each half edges to each vertex////////////////////////////////////////////
					std::pair<int, int> p_t;
					p_t.first = g[P_node_mapping[i][j]].edge[k].from;
					p_t.second = g[P_node_mapping[i][j]].edge[k].to;
					if (halfEdge_map.find(p_t)->second == i){
						Edge e1, e2;
						e1.from = j;
						e1.to = P[i].size() - 1;
						e1.elabel = elabel;
						e2.to = j;
						e2.from = P[i].size() - 1;
						e2.elabel = elabel;
						elabel--;
						/////////Here we add the half edge to the dummy node
						P[i][j].edge.push_back(e1);
						P[i][P[i].size() - 1].edge.push_back(e2);
					}
				}
			}
		}
	}
	for (unsigned int i = 0; i < P.size(); i++){
		if (P[i][P[i].size() - 1].edge.size() == 0)	P[i].pop_back();		// If a partition has no half edge at all we will delete the dummy node
	}
	/*
	for (unsigned int i = 0; i < P.size(); i++){
	for (unsigned int j = 0; j < P[i].size(); j++){
	std::cout << P[i][j].edge.size() << ' ';

	}std::cout << '\n';
	}
	std::cout << '\n';
	*/
	for (unsigned int i = 0; i < P.size(); i++){
		P[i].buildEdge();
	}
	/*
	for (int i = 0; i < b.size(); ++i) {
	std::cout << b[i] << ' ';
	}std::cout << '\n';

	
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
	cout << "P" << i << "  Sel:"<<SelectivityValue[i]<<"   ";
	for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
	cout << P_node_mapping[i][j] <<  " ";// "-" << g[P_node_mapping[i][j]].label
	}
	cout << endl;
	}
	cin.get();
*/	
	vector<int> zzzz; int f = 1;
	for (unsigned int i = 0; i < P_node_mapping.size(); i++){
		for (unsigned int j = 0; j < P_node_mapping[i].size(); j++){
//			cout << P_node_mapping[i][j] << " ";// "-" << g[P_node_mapping[i][j]].labelelse 
			
			for (int n = 0; n < zzzz.size(); n++){
				if (zzzz[n] == P_node_mapping[i][j]) { std::cout << "Wrong!!!!!!!!!" << endl; cin.get(); f = 0; }
				
			}zzzz.push_back(P_node_mapping[i][j]);
		}
		
	}
//	cout << "ZZ: ";
//	for (int cc = 0; cc < zzzz.size(); cc++){
//		std::cout << zzzz[cc] << " ";
//	}

	return P;
}
void Database::DFSCode_Inverted_Index_push(std::string &c, unsigned int graph_id, Graph &g, int m){
	std::map <std::string, Graph>::iterator pi;
	pi = Partition_inverted_index[m].find(c); 
	if (pi == Partition_inverted_index[m].end()){
		std::pair<std::string, Graph> p;
		p.first = c;
//		std::vector<int> V;
//		V.push_back(graph_id);
		p.second = g;
		p.second.graph_ids.push_back(graph_id);
		Partition_inverted_index[m].insert(p);
	}
	else{
		pi->second.graph_ids.push_back(graph_id);
	}
}

std::istream &Database::read_data(std::istream &is)
{
	int graph_counter=0;
	while (true) {
		Graph g;
		g.read(is);
		if (g.empty()) break;
		g.ID = graph_counter;
		graph_counter++;
		Graphs.push_back(g);
	}
	return is;
}

std::istream &Database::read_query(std::istream &is)
{
	while (true) {
		Graph g;
		g.read(is);
		if (g.empty()) break;
		Queries.push_back(g);
	}
	return is;
}

bool get_forward_root(Graph &g, Vertex &v, EdgeList &result)
{
	result.clear();
	for (Vertex::edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it) {
		assert(it->to >= 0 && it->to < g.size());
		if (v.label <= g[it->to].label)
			result.push_back(&(*it));
	}

	return (!result.empty());
}

Edge *get_backward(Graph &graph, Edge* e1, Edge* e2, History& history)
{
	if (e1 == e2)
		return 0;

	assert(e1->from >= 0 && e1->from < graph.size());
	assert(e1->to >= 0 && e1->to < graph.size());
	assert(e2->to >= 0 && e2->to < graph.size());

	for (Vertex::edge_iterator it = graph[e2->to].edge.begin();
		it != graph[e2->to].edge.end(); ++it)
	{
		if (history.hasEdge(it->id))
			continue;

		if ((it->to == e1->from) &&
			((e1->elabel < it->elabel) ||
			(e1->elabel == it->elabel) &&
			(graph[e1->to].label <= graph[e2->to].label)
			))
		{
			return &(*it);
		}
	}

	return 0;
}

bool get_forward_pure(Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
	result.clear();

	assert(e->to >= 0 && e->to < graph.size());

	/* Walk all edges leaving from vertex e->to.
	*/
	for (Vertex::edge_iterator it = graph[e->to].edge.begin();
		it != graph[e->to].edge.end(); ++it)
	{
		/* -e-> [e->to] -it-> [it->to]
		*/
		assert(it->to >= 0 && it->to < graph.size());
		if (minlabel > graph[it->to].label || history.hasVertex(it->to))
			continue;

		result.push_back(&(*it));
	}

	return (!result.empty());
}

bool get_forward_rmpath(Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
	result.clear();
	assert(e->to >= 0 && e->to < graph.size());
	assert(e->from >= 0 && e->from < graph.size());
	int tolabel = graph[e->to].label;

	for (Vertex::edge_iterator it = graph[e->from].edge.begin();
		it != graph[e->from].edge.end(); ++it)
	{
		int tolabel2 = graph[it->to].label;
		if (e->to == it->to || minlabel > tolabel2 || history.hasVertex(it->to))
			continue;

		if (e->elabel < it->elabel || (e->elabel == it->elabel && tolabel <= tolabel2))
			result.push_back(&(*it));
	}

	return (!result.empty());
}

void Database::IndexingMethods_push(int t){
	IndexingMethods.push_back(t);
	std::map <std::string, Graph> T_Index;
	Partition_inverted_index.push_back(T_Index);
}

void Database::Init(std::string d, std::string q, std::string o, int t, int threshold){
	ofstream oss;
	oss.open(o, fstream::app | fstream::out);
//	std::cout << "here         "  << IndexingMethods[0] << "-" << IndexingMethods[1] << "-" << IndexingMethods[2] << std::endl;
	std::ifstream isn_q(q);
	std::ifstream isn_d(d);
	read_query(isn_q);
	read_data(isn_d);
	Candidates.resize(Queries.size());
	for (unsigned int i = 0; i < Queries.size(); i++){	
		Candidates[i].resize(Graphs.size());
	}
	cout << "Loading database:" << endl;
	cout << "Number of queries:" << Queries.size() << endl;
	cout << "Number of graphs:" << Graphs.size() << endl << "Completes" << endl;

	//	For log purpose
	time_t now = time(0);
	oss << endl << ctime(&now) << "Number of Partitions:" << t << endl<<"Number of queries:" << Queries.size() << endl;
//	oss << endl << "Number of Partitions:" << t << " Threshold:" << threshold << endl << "Number of queries:" << Queries.size() << endl;
	oss << "Number of graphs:" << Graphs.size() <<  endl;

	Statistic(1);
	oss << "Avg nodes:" << TotalNodes/Graphs.size()<<" Avg edges:"<<TotalEdges/Graphs.size() << endl;

	oss.close();
}
int checkPartitionSize(std::vector < Graph > &V, int &counter,int &gap){
	if (counter > 50)
	{
		gap += 5; counter = 0;
	}
	for (int i = 0; i < V.size(); i++){
		if (V[i].size()>gap)
		{
//			std::cout << "Re-partition" << "\r";
			counter++;
			return 1;
		}
	}
	return 0;
}
void Database::Build_index(std::string d, int t, std::string o){
	ofstream oss;
	oss.open(o, fstream::app | fstream::out);
	class TSINGHUA_CLIPSE_UTIL::TimeRecorder time;
	int time1 = 0, time2 = 1;
	float t_time;
	time.check(); time1++; time2++;
	int PartitionSize = 0;
	int ProfileSize = 0;
	int GraphIDInvertedIndexSize = 0;
	
	for (unsigned int Number_Methods = 0; Number_Methods < IndexingMethods.size(); Number_Methods++){

		std::cout << "Building index" << IndexingMethods[Number_Methods] << ":" << endl;
		for (int j = 0; j < Graphs.size(); j++){
//			std::cout << j << "/" << Graphs.size() << " Size:" << Graphs[j].size() << "        Partition" << " Methods:" << Number_Methods << "\r";
			if (Graphs[j].size() < t) {
				SpecialGraphs.push_back(Graphs[j]);
				SpecialGraphsID.push_back(j);
				continue;				// If size of a graph is smaller than number of partition then we skip it
			}
			int counter = 0;
			int gap = 30;
			switch (IndexingMethods[Number_Methods]) {
			case 1:
				Partitions = Graph_partition_1(Graphs[j], t);//	std::cout << endl << "xxxxxxxxxxx  " << (IndexingMethods[Number_Methods]) << endl;

				while (checkPartitionSize(Partitions,counter,gap)){ Partitions = Graph_partition_1(Graphs[j], t); }
				break;
			case 2:
				Partitions = Graph_partition_2(Graphs[j], t);//	std::cout << endl << "xxxxxxxxxxx  " << (IndexingMethods[Number_Methods]) << endl;
				while (checkPartitionSize(Partitions, counter, gap)){ Partitions = Graph_partition_2(Graphs[j], t); }
				break;
			case 3:
				Partitions = Graph_partition_3(Graphs[j], t);//	std::cout << endl << "xxxxxxxxxxx  " << (IndexingMethods[Number_Methods]) << endl;
				while (checkPartitionSize(Partitions, counter, gap)){ Partitions = Graph_partition_3(Graphs[j], t); }
				break;
			case 4:
				Partitions = Graph_partition_4(Graphs[j], t);//	std::cout << endl << "xxxxxxxxxxx  " << (IndexingMethods[Number_Methods]) << endl;
				while (checkPartitionSize(Partitions, counter, gap)){ Partitions = Graph_partition_4(Graphs[j], t); }
				break;
			default:
				std::cout << "Unknown partition method"; std::cin.get();
			}
			
			
			for (unsigned int i = 0; i < Partitions.size(); i++){
				//std::cout << j << "/" << Graphs.size() << " Size:" << Partitions[i].size() << " Canonicalization" << "\r";
				//			Partitions[i].showGraph();
				class TSINGHUA_CLIPSE_UTIL::TimeRecorder timet;
				int a = 0, b = 1;
				float t_time;
				timet.check(); a++; b++;
	//			std::cout << endl << "Finding code.  Partition  " <<i<<"   "<< Partitions[i].size() << endl;
	//			Partitions[i].showGraph();
				PSS::DFSCode c = CanonicalCode(Partitions[i]); 
				timet.check();
				double slow = timet.diffTime(a++, b++);
				if (slow>20){
					Partitions[i].printGraph("workspace/output/Slowgraph.debug");
	//				std::cout <<slow<< " ";
				}
				string t = c.toString();
				Partitions[i].buildProfile();
				DFSCode_Inverted_Index_push(t, j, Partitions[i], Number_Methods); 
				PartitionSize += t.size();
				ProfileSize += (Partitions[i].Label_statistic_v.size() + Partitions[i].Label_statistic_e.size()) * 8;				
			}
			Partitions.resize(0);
		}
		time.check();
		t_time = time.diffTime(time1++, time2++);
		std::cout << "Time for building index:" << t_time << ". Index size:" << Partition_inverted_index[Number_Methods].size() << endl;

		for (std::map <std::string, Graph>::iterator pi = Partition_inverted_index[Number_Methods].begin(); pi != Partition_inverted_index[Number_Methods].end(); pi++){
			GraphIDInvertedIndexSize += (pi->second.graph_ids.size() * 4);
		}

		switch (IndexingMethods[Number_Methods]) {
		case 1:
			oss << "Time for building index (BFS):" << t_time << ". Index size:" << Partition_inverted_index[Number_Methods].size() << " Space:" << (double)(PartitionSize + ProfileSize+GraphIDInvertedIndexSize) / (double)1048576 << "MB" << endl;
			break;
		case 2:
			oss << "Time for building index (DFS):" << t_time << ". Index size:" << Partition_inverted_index[Number_Methods].size() << " Space:" << (double)(PartitionSize + ProfileSize + GraphIDInvertedIndexSize) / (double)1048576 << "MB" << endl;
			break;
		case 3:
			oss << "Time for building index (Selectivity1):" << t_time << ". Index size:" << Partition_inverted_index[Number_Methods].size() << " Space:" << (double)(PartitionSize + ProfileSize + GraphIDInvertedIndexSize) / (double)1048576 << "MB" << endl;
			break;
		case 4:
			oss << "Time for building index (Selectivity2):" << t_time << ". Index size:" << Partition_inverted_index[Number_Methods].size() << " Space:" << (double)(PartitionSize + ProfileSize + GraphIDInvertedIndexSize) / (double)1048576 << "MB" << endl;
			break;
		default:
			std::cout << "Unknown partition method"; std::cin.get();
		}
	}
	oss << "Number of special graphs:" << SpecialGraphsID.size() << endl << endl;;
	oss.close();
}

void candidates_updates(std::vector<std::vector<int>> &result, std::vector<std::vector<int>> &Candidates){
	for (unsigned int i = 0; i < Candidates.size(); i++){
		for (unsigned int j = 0; j < Candidates[i].size(); j++){
			if ((Candidates[i][j] == 0) && ((result[i][j] == 1)))Candidates[i][j] = 1;
			if (result[i][j] == -1)Candidates[i][j] = -1;
		}
	}
}

void Database::RemoveGraph(int size){
	vector<Graph> temp;
	for (int i = 0; i < Graphs.size(); i++){
		if (Graphs[i].size() <= size){ temp.push_back(Graphs[i]); }
	}
	Graphs.resize(0);
	Graphs = temp;
}
void  Database::Random_makeup_queries(int n){
	Queries.resize(0);
	std::vector<int> a(Graphs.size());
	for (int i = 0; i<Graphs.size(); ++i) a[i] = i;
	//	srand(time(NULL));
	random_unique(a.begin(), a.end(), n);
	//	random_unique2(a.begin(), a.end(), a.size());
	/*
	for (int i = 0; i < a.size(); ++i) {
	std::cout << a[i] << ' ';
	}std::cout << '\n';
	*/
	
	for (int i = 0; i<n; ) {
		if (Graphs[a[i]].size() < 101){
			Graphs[a[i]].printGraph("workspace/input/100SampleProtein2");
			i++;
		}
		
	}
}

void Database::Statistic(int show){
	
	for (int i = 0; i < Graphs.size(); i++){
		for (int j = 0; j < Graphs[i].size(); j++){
			TotalNodes++;
			std::map <int,int>::iterator l;
			int c = Graphs[i][j].label;
			l = Label_statistic_v.find(c);
			if (l == Label_statistic_v.end()){
				std::pair<int,int> p;
				p.first = c;
				p.second = 1;
				Label_statistic_v.insert(p);
			}
			else{
				l->second++;
			}
			for (int k = 0; k < Graphs[i][j].edge.size(); k++){
				TotalEdges++;
				std::map <int, int>::iterator l;
				int c = Graphs[i][j].edge[k].elabel;
				l = Label_statistic_e.find(c);
				if (l == Label_statistic_e.end()){
					std::pair<int, int> p;
					p.first = c;
					p.second = 1;
					Label_statistic_e.insert(p);
				}
				else{
					l->second++;
				}
			}
		}
	}

	if (show){
		std::map <int, int>::iterator l;
		for (l = Label_statistic_v.begin(); l != Label_statistic_v.end(); l++){
			std::cout << "V_Label:" << l->first << " freqency:" << l->second << endl;
		}
		for (l = Label_statistic_e.begin(); l != Label_statistic_e.end(); l++){
			std::cout << "E_Label:" << l->first << " freqency:" << l->second << endl;
		}
		std::cout << "Avg nodes:" << TotalNodes / Graphs.size() << " Avg edges:" << TotalEdges / Graphs.size() << endl;
	}

}
int profilePrune(Graph &a, Graph &b){
	std::map <int,int>::iterator p1,p2;
	p1 = a.Label_statistic_v.begin();
	p2 = b.Label_statistic_v.begin();
	while ((p1 != a.Label_statistic_v.end()) && (p2 != b.Label_statistic_v.end())){
		if (p1->first == -2) {
			p1++; continue;
		}
		if (p1->first > p2->first) p2++;
		else if (p1->first < p2->first) return 0;
		else if (p1->second > p2->second) return 0;
		else { p1++; p2++; }
	}
	return 1;
}
void Database::Query_processing(int t, int threshold, std::string o,int profile){
	std::cout << "Processing query:" << endl;
	int ProfilePruned = 0;
	int NumberofSubiso = 0;
//	std::vector<unsigned int> result;
	std::map <std::string, Graph>::iterator pi;
	ofstream oss;
	oss.open(o, fstream::app | fstream::out);
	class TSINGHUA_CLIPSE_UTIL::TimeRecorder time;
	int time1 = 0, time2 = 1;
	float t_time;
	time.check(); time1++; time2++;
	for (unsigned int i = 0; i < Queries.size(); i++){
		Queries[i].buildProfile();
	}
	///////////////////	// For testing best partition number
	vector<int> avgCandidate(Queries.size(),0);	
	Candidates.resize(Queries.size());
	for (unsigned int i = 0; i < Queries.size(); i++){
		Candidates[i].resize(Graphs.size());
	}
	/////////////////	// For testing best partition number
	for (unsigned int Number_Methods = 0; Number_Methods < IndexingMethods.size(); Number_Methods++){
		std::vector<std::vector<int>> result;
		result.resize(Queries.size());
		for (unsigned int i = 0; i < Queries.size(); i++){
			result[i].resize(Graphs.size());
		}
		for (unsigned int i = 0; i < (*this).Queries.size(); i++){
			int j = 0;
			int xx = 0;
			int k = 0;
			for (pi = Partition_inverted_index[Number_Methods].begin(); pi != Partition_inverted_index[Number_Methods].end(); pi++){
				//		std::cout << "Query " << i << ": " << k << "/" << Partition_inverted_index[Number_Methods].size() << " Methods:" << Number_Methods << "\r"; k++;
				if (profile){
					NumberofSubiso++;
					if (profilePrune(pi->second, Queries[i]))
					{
						ull_ge::halfedge_ulliso u(pi->second, Queries[i]);
						if (u.match()){
							for (unsigned int j = 0; j < pi->second.graph_ids.size(); j++){
								result[i][pi->second.graph_ids[j]]++;
							}
						}
					}
					else{
						ProfilePruned++;
					}
				}
				else
				{
					ull_ge::halfedge_ulliso u(pi->second, Queries[i]);
					if (u.match()){
						for (unsigned int j = 0; j < pi->second.graph_ids.size(); j++){
							result[i][pi->second.graph_ids[j]]++;
						}
					}
				}
			}
		}

		if (0){		// Debug
			for (int k = 0; k < result.size(); k++){
				for (int l = 0; l < result[k].size(); l++){
					std::cout << result[k][l] << " ";
					if (result[k][l]>(t - threshold))result[k][l] = 1;
					else result[k][l] = 0;
				}std::cout << std::endl;
			}
			std::cout << "------------------------------" << std::endl;
			for (int k = 0; k < result.size(); k++){
				for (int l = 0; l < result[k].size(); l++){
					std::cout << result[k][l] << " ";
					if (result[k][l] == 1){ oss << l << " "; }
				}std::cout << std::endl;
			}
		}

		for (int k = 0; k < result.size(); k++){
			for (int l = 0; l < result[k].size(); l++){
				if (result[k][l]>(t - threshold))result[k][l] = 1;
				else result[k][l] = -1;
			}
		}

		candidates_updates(result, Candidates);
	}
	for (int k = 0; k < Candidates.size(); k++){
		oss << "Query " << k << "(size"<<Queries[k].size()<<")'s Candidates:" ;
		int i = 0;
		for (int l = 0; l < Candidates[k].size(); l++){
			if (Candidates[k][l] == 1){ i++; }
		}
		oss << " (Pruning ability:" << i << "/" << Candidates[k].size() << "=" << (double)i / (double)Candidates[k].size() << ")" << endl;
//		oss << "Pool " << k << " :";
		avgCandidate[k] = i;
//		for (int l = 0; l < Candidates[k].size(); l++){
//			if (Candidates[k][l] == 1){ oss << l << " "; }
//		}
//		oss << endl;
	}
	time.check(); 
	t_time = time.diffTime(time1++, time2++);
	///////////////////	// For testing best partition number
	double CanSize = 0;
	for (unsigned int i = 0; i < avgCandidate.size(); i++){
		CanSize+=avgCandidate[i];
	}
	CanSize = (double)CanSize / (double)Queries.size();
	oss << "Threshold: "<<threshold<<"  Number of partitions: "<<t<< endl;
	oss << "(Average Pruning ability:" << CanSize << "/" << Graphs.size() << "=" << (double)CanSize / (double)Graphs.size() << ")" << endl;
	///////////////////	// For testing best partition number

	oss << "Total query processing time:" << t_time << endl;
	if (profile==1)
		oss << "Profile pruned Iso computation Ratio:" << ProfilePruned << "/" << NumberofSubiso<<": "<<(double)ProfilePruned / (double)NumberofSubiso << endl;
	else
	oss << "Graph profile not used" << endl;
	cout << "Total query processing time:" << t_time << endl;
	
	GED::CostFunction cf;
	GED::EditDistance ed;
	for (int k = 0; k < Candidates.size(); k++){
		//oss << "Query " << k << "'s verfied result:" << endl;
		int falsepositive = 0;int can = 0;
		for (int l = 0; l < Candidates[k].size(); l++){
			if (Candidates[k][l] == 1)
			{
				//oss << l << '\t';
				can++;
			}
		}
		//oss << endl;
		for (int l = 0; l < Candidates[k].size(); l++){
			int realGED = ed.getEditDistance(Queries[k], Graphs[l], cf, threshold + 1, 0);
			if (realGED>threshold) falsepositive++;
			//oss << realGED << '\t';
		}
		//oss << endl;
		oss << "False positive rate:" << falsepositive << "/" << can << "=" << (double)falsepositive / (double)can << endl;
	}
	time.check(); 
	t_time = time.diffTime(time1++, time2++);
	oss << "Total query verification time:" << t_time << endl;
	cout << "Total query verification time:" << t_time << endl;
	cout << "Results are saved to:" << o << endl;
	oss.close();
}

DFSCode Database::CanonicalCode(Graph &GRAPH_IS_MIN){
	DFSCode DFS_CODE_IS_MIN;

	if (GRAPH_IS_MIN.size() == 1)
	{
		DFSCode dummy;
		DFS t;
		t.from = 0;
		t.to = 0;
		t.fromlabel = GRAPH_IS_MIN[0].label;
		t.tolabel = GRAPH_IS_MIN[0].label;
		t.elabel = -1;
		dummy.push_back(t);
		return dummy;
	}

	Projected_map3	root;
	EdgeList        edges;

	for (unsigned int from = 0; from < GRAPH_IS_MIN.size(); ++from)
	if (get_forward_root(GRAPH_IS_MIN, GRAPH_IS_MIN[from], edges))
	for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
	{
		root[GRAPH_IS_MIN[from].label][(*it)->elabel][GRAPH_IS_MIN[(*it)->to].label].push(0, *it, 0);
//		std::cout << (*it)->from << "-" << (*it)->to << std::endl;
	}

	///////////Debug/////////////////
	if (0)
	{
		using namespace std;
		cout <<endl<< "------Debug mincode------------ " << endl;
		for (Projected_iterator3 i=root.begin(); i != root.end(); i++){
			cout << "root from label " << i->first <<endl;
			for (Projected_iterator2 j = i->second.begin(); j != i->second.end(); j++){
				cout << "	root e label " << j->first << endl;
				for (Projected_iterator1 k = j->second.begin(); k != j->second.end(); k++){
					cout << "		root to label " << k->first << endl;
				}
			}
		}
		cout << "-------------------------------------- " << endl;
	}
	///////////Debug////////////////

	Projected_iterator3 fromlabel = root.begin();
	Projected_iterator2 elabel = fromlabel->second.begin();
	Projected_iterator1 tolabel = elabel->second.begin();

	DFS_CODE_IS_MIN.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
	
	return project_is_min(tolabel->second, DFS_CODE_IS_MIN, GRAPH_IS_MIN);
}



DFSCode Database::project_is_min(Projected &projected, DFSCode &DFS_CODE_IS_MIN, Graph &GRAPH_IS_MIN)
{
	const RMPath& rmpath = DFS_CODE_IS_MIN.buildRMPath();

	if (0)
	{	
		using namespace std;
		std::cout << endl << "Debug-----rmpath-------------------------------" << endl;
		for (unsigned int i = 0; i < rmpath.size(); i++){
		std::cout << rmpath[i] << " ";		
		}
		cout << endl << "----------------------------------------------------" << endl;
	}

	int minlabel = DFS_CODE_IS_MIN[0].fromlabel;
	int maxtoc = DFS_CODE_IS_MIN[rmpath[0]].to;

	if (0)
	{
		using namespace std;
		std::cout << endl << "minlabel:" << minlabel <<" maxtoc:"<<maxtoc<< endl;
	}

	{
		Projected_map1 root;
		bool flg = false;
		int newto = 0;

		for (int i = rmpath.size() - 1; !flg && i >= 1; --i) {
			for (unsigned int n = 0; n < projected.size(); ++n) {
				PDFS *cur = &projected[n];
				History history(GRAPH_IS_MIN, cur);
				Edge *e = get_backward(GRAPH_IS_MIN, history[rmpath[i]], history[rmpath[0]], history);
				if (e) {
					if (0)
					{//	if there is backward edge
						std::cout << "backward e:" << e->from << "-" << e->to << std::endl;
					}
					
					root[e->elabel].push(0, e, cur);
					newto = DFS_CODE_IS_MIN[rmpath[i]].from;
					flg = true;
				}
			}
		}

		if (flg) {
			Projected_iterator1 elabel = root.begin();
			DFS_CODE_IS_MIN.push(maxtoc, newto, -1, elabel->first, -1);
//			if (DFS_CODE[DFS_CODE_IS_MIN.size() - 1] != DFS_CODE_IS_MIN[DFS_CODE_IS_MIN.size() - 1]) return false;
			if (0)
			{
				DFS_CODE_IS_MIN.showDFScode();
			}
			return project_is_min(elabel->second, DFS_CODE_IS_MIN, GRAPH_IS_MIN);
		}
	}

	{
		bool flg = false;
		int newfrom = 0;
		Projected_map2 root;
		EdgeList edges;
		
//		std::cout << " projected.size() " << projected.size() << std::endl;

		for (unsigned int n = 0; n < projected.size(); ++n) {
			PDFS *cur = &projected[n];
//			std::cout << " *cur id:" << cur->id << std::endl; cur->edge->show();
			History history(GRAPH_IS_MIN, cur);
			if (get_forward_pure(GRAPH_IS_MIN, history[rmpath[0]], minlabel, history, edges)) {
				flg = true;
				newfrom = maxtoc;
				for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
					root[(*it)->elabel][GRAPH_IS_MIN[(*it)->to].label].push(0, *it, cur);
			}
		}

		for (int i = 0; !flg && i < (int)rmpath.size(); ++i) {
			for (unsigned int n = 0; n < projected.size(); ++n) {
				PDFS *cur = &projected[n];
				History history(GRAPH_IS_MIN, cur);
				if (get_forward_rmpath(GRAPH_IS_MIN, history[rmpath[i]], minlabel, history, edges)) {
					flg = true;
					newfrom = DFS_CODE_IS_MIN[rmpath[i]].from;
					for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
						root[(*it)->elabel][GRAPH_IS_MIN[(*it)->to].label].push(0, *it, cur);
				}
			}
		}

		if (flg) {
			Projected_iterator2 elabel = root.begin();
			Projected_iterator1 tolabel = elabel->second.begin();
			DFS_CODE_IS_MIN.push(newfrom, maxtoc + 1, -1, elabel->first, tolabel->first);
//			if (DFS_CODE[DFS_CODE_IS_MIN.size() - 1] != DFS_CODE_IS_MIN[DFS_CODE_IS_MIN.size() - 1]) return false;
			if (0)
			{
				DFS_CODE_IS_MIN.showDFScode();
			}
			return project_is_min(tolabel->second, DFS_CODE_IS_MIN, GRAPH_IS_MIN);
		}
	}

	return DFS_CODE_IS_MIN;
}

}
