#ifndef Database_h
#define Database_h

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <string>
#include <strstream>
#include <unordered_map>
#include "Matrix.h"
namespace PSS{

class Edge{
public:
	int from;
	int to;
	int elabel;
	unsigned int id;
	void show(){
		std::cout << "edge:" << id << " " << from << "-" << elabel << "-" << to << std::endl;
	}
};

class Vertex{
public:
	int v_id;
	int label;
	std::vector<Edge> edge;
	typedef std::vector<Edge>::iterator edge_iterator;
	void push(int from,int to,int elabel){
		edge.resize(edge.size() + 1);
		edge[edge.size() - 1].from = from;
		edge[edge.size() - 1].to = to;
		edge[edge.size() - 1].elabel = elabel;
		return;
	}	
	void showEdge(){
		std::cout << std::endl;
		for (unsigned int i = 0; i < edge.size(); i++){
			std::cout << edge[i].from <<"-"<< edge[i].to << std::endl;
		}
		std::cout << std::endl;
	}
};

class Graph : public std::vector<Vertex>{
private:
	unsigned int edge_size_;
	
public:
//	vector <int> Subgraph_Vertex_mapping;	// This is to map the ID of nodes in a subgraph to the Id of nodes in the subgraph's supergraph 
	typedef std::vector<Vertex>::iterator vertex_iterator;
	int ID=-1;									// The ID of this graph in the database
	std::map<int, int>													 Label_statistic_v;
	std::map<int, int>													 Label_statistic_e;
	bool directed;
	unsigned int vertex_size() { return (unsigned int)size(); }
	unsigned int edge_size()   { return edge_size_; }
	std::vector<int> graph_ids;		// For partition inverted index
	Graph(){
		directed = false; edge_size_ = 0;
	}
	std::istream &read(std::istream &); // read
	int buildEdge();
	void buildProfile();				//
	void showGraph();
	bool hasEdge(int i,int j){
		for (unsigned int k = 0; k < (*this)[i].edge.size(); k++){
			if ((*this)[i].edge[k].to == j) return 1;
		}
		return 0;
	}
	bool hasEdgeElabel(int i, int j, int ela){
		for (unsigned int k = 0; k < (*this)[i].edge.size(); k++){
			if (((*this)[i].edge[k].to == j) && ((*this)[i].edge[k].elabel == ela)) return 1;
		}
		return 0;
	}
	Matrix<int> ToAdjacencyMatrix(){
		Matrix<int> adjacencyMatrix;
		adjacencyMatrix.resize(this->size(),this->size());
		for (unsigned int i = 0; i < this->size(); i++){
			for (unsigned int j = 0; j < (*this)[i].edge.size(); j++){
				adjacencyMatrix.m_matrix[(*this)[i].edge[j].from][(*this)[i].edge[j].to] = 1;
			}
		}
		return adjacencyMatrix;
	}
	void Assign_vertexid(){
		for (unsigned int i = 0; i < this->size(); i++) (*this)[i].v_id = i;
	}
	Graph remove_vertex(int index){		//Remove a node in a graph and all the edges connected to this node, then return the new graph, subgraph_vertex_mapping record the mapping of nodes in these two graphs
		Graph new_graph;

		unsigned int j=0;
		for (unsigned int i = 0; i < this->size(); i++)
		{
			Vertex v = (*this)[i];
			v.edge.resize(0);
			if (i == index) continue;
			else {
//				new_graph.Subgraph_Vertex_mapping.push_back((*this).Subgraph_Vertex_mapping[i]);
				j++;
				for (unsigned int k = 0; k < (*this)[i].edge.size(); k++){
					Edge e = (*this)[i].edge[k];
					if (e.to = index) { continue; }
					if (e.to > index)	{ e.to = e.to - 1; }
					if (e.from > index)	{ e.from = e.from - 1; }
					v.push(e.from, e.to,e.elabel);
				}
				new_graph.push_back(v);
			}			
		}
		new_graph.buildEdge();
		return new_graph;
	}

	void printGraph(std::string d);
	void printGraph(std::string d, std::string code);
};


class DFS {
public:
	int from;
	int to;
	int fromlabel;
	int elabel;
	int tolabel;
	friend bool operator == (const DFS &d1, const DFS &d2)
	{
		return (d1.from == d2.from && d1.to == d2.to && d1.fromlabel == d2.fromlabel
			&& d1.elabel == d2.elabel && d1.tolabel == d2.tolabel);
	}
	friend bool operator != (const DFS &d1, const DFS &d2) { return (!(d1 == d2)); }
	DFS() : from(0), to(0), fromlabel(0), elabel(0), tolabel(0) {};
};

typedef std::vector<int> RMPath;

class DFSCode : public std::vector<DFS>{
private:
	RMPath rmpath;
public:
	const RMPath& buildRMPath()
	{
		rmpath.clear();

		int old_from = -1;

		for (int i = size() - 1; i >= 0; --i) {
//			std::cout <<"buildRMPath "<< (*this)[i].from << std::endl;
			if ((*this)[i].from < (*this)[i].to && // forward
				(rmpath.empty() || old_from == (*this)[i].to))
			{				
				rmpath.push_back(i);
				old_from = (*this)[i].from;
			}
		}

		return rmpath;
	}

	void push(int from, int to, int fromlabel, int elabel, int tolabel)
	{
		resize(size() + 1);
		DFS &d = (*this)[size() - 1];

		d.from = from;
		d.to = to;
		d.fromlabel = fromlabel;
		d.elabel = elabel;
		d.tolabel = tolabel;
	}
	void pop() { resize(size() - 1); }

	std::string toString(){
		std::string code;
		
		std::vector<char> code_char;
		for (int i = 0; i < (*this).size(); i++){
			code_char.push_back((*this)[i].from);
			code_char.push_back((*this)[i].to);
			code_char.push_back((*this)[i].fromlabel);
			code_char.push_back((*this)[i].elabel);
			code_char.push_back((*this)[i].tolabel);
		}
		for (unsigned int i = 0; i < code_char.size(); i++){
//			code << code_char[i];
			code += std::to_string(code_char[i]);
		}
//		std::cout << code << std::endl;
		return code;
	}

	bool toGraph(PSS::Graph &g){
		g.clear();

		for (DFSCode::iterator it = begin(); it != end(); ++it) {
			g.resize(std::max(it->from, it->to) + 1);

			if (it->fromlabel != -1)
				g[it->from].label = it->fromlabel;
			if (it->tolabel != -1)
				g[it->to].label = it->tolabel;

			g[it->from].push(it->from, it->to, it->elabel);
			if (g.directed == false)
				g[it->to].push(it->to, it->from, it->elabel);
		}
		g.buildEdge();
		return (true);
	}
	PSS::Graph  toGraph2(){
		PSS::Graph g;

		for (DFSCode::iterator it = begin(); it != end(); ++it) {
			g.resize(std::max(it->from, it->to) + 1);

			if (it->fromlabel != -1)
				g[it->from].label = it->fromlabel;
			if (it->tolabel != -1)
				g[it->to].label = it->tolabel;

			g[it->from].push(it->from, it->to, it->elabel);
			if (g.directed == false)
				g[it->to].push(it->to, it->from, it->elabel);
		}
		g.buildEdge();
		return g;
	}
	void showDFScode(){
		using namespace std;
		cout << "-----------------------------" << endl;
//		std::cout << endl << "-----DFScode-------------------" << endl;
		if ((*this).size() == 0)	std::cout << "A dummy DFScode. (Graph only has one node)" << std::endl;
		for (unsigned int i = 0; i < size(); i++){
			std::cout << (*this)[i].from << " " << (*this)[i].to << " " << (*this)[i].fromlabel << " " << (*this)[i].elabel << " " << (*this)[i].tolabel << std::endl;
		}
//		cout  << "-----------------------------" << endl;
	}
};





struct PDFS {
	unsigned int id;	// ID of the original input graph
	Edge        *edge;
	PDFS        *prev;
	PDFS() : id(0), edge(0), prev(0) {};
};

class History : public std::vector<Edge*> {
private:
	std::vector<int> edge;
	std::vector<int> vertex;

public:
	bool hasEdge(unsigned int id) { return (bool)edge[id]; }
	bool hasVertex(unsigned int id) { return (bool)vertex[id]; }
	void History::build(Graph &graph, PDFS *e)
	{
		// first build history
		clear();
		edge.clear();
		edge.resize(graph.edge_size());
		vertex.clear();
		vertex.resize(graph.size());

		if (e) {
			push_back(e->edge);
			edge[e->edge->id] = vertex[e->edge->from] = vertex[e->edge->to] = 1;

			for (PDFS *p = e->prev; p; p = p->prev) {
				push_back(p->edge);	// this line eats 8% of overall instructions(!)
				edge[p->edge->id] = vertex[p->edge->from] = vertex[p->edge->to] = 1;
			}
			std::reverse(begin(), end());
		}
	}
	History() {};
	History(Graph& g, PDFS *p) { build(g, p); }

};

class Projected : public std::vector<PDFS> {
public:
	void push(int id, Edge *edge, PDFS *prev)
	{
		resize(size() + 1);
		PDFS &d = (*this)[size() - 1];
		d.id = id; d.edge = edge; d.prev = prev;
	}
};

typedef std::vector <Edge*> EdgeList;

bool  get_forward_pure(Graph&, Edge *, int, History&, EdgeList &);
bool  get_forward_rmpath(Graph&, Edge *, int, History&, EdgeList &);
bool  get_forward_root(Graph&, Vertex&, EdgeList &);
Edge *get_backward(Graph&, Edge *, Edge *, History&);

class Database{
private:
	typedef std::map<int, std::map <int, std::map <int, Projected> > >           Projected_map3;
	typedef std::map<int, std::map <int, Projected> >                            Projected_map2;
	typedef std::map<int, Projected>                                             Projected_map1;
	typedef std::map<int, std::map <int, std::map <int, Projected> > >::iterator Projected_iterator3;
	typedef std::map<int, std::map <int, Projected> >::iterator                  Projected_iterator2;
	typedef std::map<int, Projected>::iterator                                   Projected_iterator1;
	typedef std::map<int, std::map <int, std::map <int, Projected> > >::reverse_iterator Projected_riterator3;

	std::map<int, int>													 Label_statistic_v;
	std::map<int, int>													 Label_statistic_e;
	int TotalNodes = 0;
	int TotalEdges = 0;
	std::istream &read_data(std::istream &);
	std::istream &read_query(std::istream &);
	std::vector<int> IndexingMethods;

public:
	std::vector < Graph >       Graphs;
	std::vector < Graph >       SpecialGraphs;		// This is used to store graphs that can not be processed by the frame work for example, we can not partition a 2 nodes graph to three partitions
	std::vector < int >			SpecialGraphsID;
	std::vector < Graph >       Queries;
	std::vector<DFSCode>		DFSCodes;
	std::vector < Graph >       Partitions;			// may have half edge
	std::vector<DFSCode>		Partition_DFSCodes; // may have half edge
	std::vector<std::map <std::string, Graph>>			Partition_inverted_index;
	std::vector<std::vector<int>> Candidates;
	
	void IndexingMethods_push(int t);
	void DFSCode_Inverted_Index_push(std::string &c, unsigned int graph_id, Graph &, int m);
	DFSCode minCode;
	DFSCode CanonicalCode(Graph &GRAPH_IS_MIN);
	DFSCode project_is_min(Projected &projected, DFSCode &DFS_CODE_IS_MIN, Graph &GRAPH_IS_MIN);
	void Random_makeup_queries(int n);
	std::vector < Graph > Graph_partition_1(Graph &, unsigned int);			//Random Partion BFS
	std::vector < Graph > Graph_partition_2(Graph &, unsigned int);			//Random Partion DFS
	std::vector < Graph > Graph_partition_3(Graph &, unsigned int);			//Selectivity Size Partition1 divide
	std::vector < Graph > Graph_partition_4(Graph &, unsigned int);			//Selectivity Size Partition2 product
	void Build_index(std::string d, int t, std::string o);								//Building index for all methods
	void Init(std::string d, std::string q, std::string o, int t, int threshold);
	void Query_processing(int t, int threshold, std::string o,int profile);			// Random partition
	void Statistic(int show);
	void RemoveGraph(int size);
};






}

#endif