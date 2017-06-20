#ifndef TreeNode_h
#define TreeNode_h
#include <iostream>
#include <vector>
#include "Database.h"
#include <map>
#include <string>
#include <vector>
#include <set>
#include <algorithm> 
#include "Hungarian.h"
#include "GraphEditDistance.h"
namespace GED{
class CostFunction{
private: 	
	/**
	* the constant cost for node and edge deletions/insertions
	*/
	double nodeCost;
	double edgeCost;
public:
	CostFunction(){ nodeCost = 1; edgeCost = 1; }
	CostFunction(double n, double e){ nodeCost = n; edgeCost =e; }
	/**
	* @return the constant cost for node/edge substitution
	*/
	double getCost(PSS::Vertex u, PSS::Vertex v) {
		double cost = 0;
		if (u.label == v.label){return 0;}
		else { return 1; }
	}
	double getCost(PSS::Edge u, PSS::Edge v) {
		double cost = 0;
		if (u.elabel == v.elabel){ return 0; }
		else { return 1; }
	}
	/**
	* @return the constant cost for node deletion/insertion
	*/
	double getNodeCosts() {
		return nodeCost;
	}
	double getEdgeCosts() {
		return edgeCost;
	}
};


class MatrixGenerator {
private:
	/**
	* the cource and target graph whereon the cost matrix is built
	*/
	PSS::Graph source, target;
	/**
	* the cost function actually employed
	*/
	CostFunction cf;
	std::string adj;
	std::vector<PSS::Edge> edges1;
	std::vector<PSS::Edge> edges2;
public:
	MatrixGenerator() {
		adj = "best";
	}
	/**
	* constructs a MatrixGenerator
	* @param costFunction
	* @param outputCostMatrix
	*/
	MatrixGenerator(CostFunction costFunction) {
		this->cf = costFunction;
		adj = "worst";
	}
	/**
	* @param centrality
	* @return the cost matrix for two graphs @param sourceGraph and @param targetGraph
	* |         |
	* | c_i,j   | del
	* |_________|______
	* |         |
	* |  ins    |	0
	* |         |
	*
	*/
	Matrix<double> getMatrix(PSS::Graph sourceGraph, PSS::Graph targetGraph) {
		this->source = sourceGraph;
		this->target = targetGraph;
		int sSize = sourceGraph.size();
		int tSize = targetGraph.size();
		int dim = sSize + tSize;
		Matrix<double> matrix(dim, dim), edgeMatrix;

		PSS::Vertex u, v;
		for (unsigned int i = 0; i < sSize; i++) {		//upper left corner
			u = source[i];
			for (unsigned int j = 0; j < tSize; j++) {
				v = target[j];
				double costs = cf.getCost(u, v);
				if (adj == "worst"){
					costs += u.edge.size()*cf.getEdgeCosts();
					costs += v.edge.size()*cf.getEdgeCosts();
				}
				if (adj == "best"){
					Matrix<double> edgeMatrix = this->getEdgeMatrix(u, v);
					Hungarian ha(edgeMatrix);
					ha.max_to_min();
					costs += ha.solve();

				}
				matrix.m_matrix[i][j] = costs;
			}
		}

		for (unsigned int i = sSize; i < dim; i++) {		//upper right corner
			for (unsigned int j = 0; j < tSize; j++) {
				if ((i - sSize) == j) {
					v = target[j];
					double costs = cf.getNodeCosts();
					if ((adj == "worst") || (adj == "best")){
						double f = v.edge.size();
						costs += (f * cf.getEdgeCosts());
						matrix.m_matrix[i][j] = costs;
					}
				}
				else {
					matrix.m_matrix[i][j] = 1000;
				}
			}
		}

		for (unsigned int i = 0; i < sSize; i++) {		// lower left
			u = source[i];
			for (unsigned int j = tSize; j < dim; j++) {
				if ((j - tSize) == i) {
					double costs = cf.getNodeCosts();;
					if ((adj == "worst") || (adj == "best")){
						double f = u.edge.size();
						costs += (f * cf.getEdgeCosts());
					}
					matrix.m_matrix[i][j] = costs;
				}
				else {
					matrix.m_matrix[i][j] = 1000;
				}
			}
		}

		for (unsigned int i = sSize; i < dim; i++) {		// lower right
			for (unsigned int j = tSize; j < dim; j++) {
				matrix.m_matrix[i][j] = 0.0;
			}
		}
		return matrix;
	}

	Matrix<double> getEdgeMatrix(PSS::Vertex u, PSS::Vertex v) {
		int uSize = u.edge.size();
		int vSize = v.edge.size();
		int dim = uSize + vSize;
		Matrix<double> edgeMatrix(dim, dim);
		PSS::Edge e_u;
		PSS::Edge e_v;
		for (int i = 0; i < uSize; i++) {
			e_u = u.edge[i];
			for (int j = 0; j < vSize; j++) {
				e_v = v.edge[j];
				double costs = cf.getCost(e_u, e_v);
				edgeMatrix.m_matrix[i][j] = costs;
			}
		}
		for (int i = uSize; i < dim; i++) {
			for (int j = 0; j < vSize; j++) {
				// diagonal
				if ((i - uSize) == j) {
					double costs = cf.getEdgeCosts();
					edgeMatrix.m_matrix[i][j] = costs;
				}
				else {
					edgeMatrix.m_matrix[i][j] = 1000;
				}
			}
		}
		for (int i = 0; i < uSize; i++) {
			for (int j = vSize; j < dim; j++) {
				// diagonal
				if ((j - vSize) == i) {
					double costs = cf.getEdgeCosts();
					edgeMatrix.m_matrix[i][j] = costs;
				}
				else {
					edgeMatrix.m_matrix[i][j] = 1000;
				}
			}
		}
		return edgeMatrix;
	}
};

class TreeNode{
private:
	std::vector<int> matching;						/** nodes of g1 are mapped to...*/
	std::vector<int> inverseMatching;				/** nodes of g2 are mapped to...*/
	std::vector<std::vector<unsigned int>> doubleMatching;
	double cost;											/** the current cost of this partial solution*/
	CostFunction cf;										/** the cost function defines the cost of individual node/edge operations*/
	double factor;											/** weighting factor for edge operations = 0.5 if undirected edges are used (1.0 otherwise)*/

	std::vector<double> costMatrix;
	PSS::Graph originalGraph1;								/** the original graphs */
	PSS::Graph originalGraph2;
	int depth = 0;
	PSS::Graph unusedNodes1;								/** the graphs where the processed nodes are removed */
	PSS::Graph unusedNodes2;
//	EditDistance editDistance;
	/**
	* @return number of adjacent edges of node with index @param i
	* NOTE: only edges (i,j) are counted if
	* j-th node hae been processed (deleted or substituted)
	*/
	int getNumberOfAdjacentEdges(std::vector<int> m, PSS::Graph a, int i) {	// This function is for undirected graph only, we only count edges that from i to j, not count from j to i
		int e = 0;
		for (int j = 0; j < a[i].edge.size(); j++){
			if (m[a[i].edge[j].to] != -1){ // count edges only if other end has been processed
					e += 1;
			}
		}
		return e;
	}

	void addCost(double c) {
		cost += c;
	}

	void processEdges(TreeNode &tn, PSS::Vertex start, int startIndex, PSS::Vertex end, int endIndex) {
		
		std::vector<int> unconnectedNode(originalGraph1.size(), 0);		// We set a node to 0 if it is not connected to startIndex
		unconnectedNode[startIndex] = 1; 
		for (int e = 0; e < tn.originalGraph1[startIndex].edge.size(); e++){// there is an edge between start and start2
			PSS::Edge edge = tn.originalGraph1[startIndex].edge[e]; 
			int start2Index = edge.to;
			unconnectedNode[start2Index] = 1;		//We set nodes to 1 if they are connected to startIndex
			if (tn.matching[start2Index] != -1) { // other end has been handled
				int end2Index = tn.matching[start2Index]; 
		//		std::cout << " processEdges s-s':" <<startIndex<<  "-" <<start2Index<< " t-t':" <<  endIndex << "-" << end2Index << " size of edges of s " << tn.originalGraph1[startIndex].edge.size() << std::endl;
				if (end2Index >= 0) {		//substitution
					if (tn.originalGraph2.hasEdge(endIndex,end2Index)) {	// if there is an edge in graph2, we deal with the edge substitution
						PSS::Edge edge2;
						for (unsigned int i= 0; i < tn.originalGraph2[endIndex].edge.size(); i++){
							if (end2Index == originalGraph2[endIndex].edge[i].to)
								edge2 = tn.originalGraph2[endIndex].edge[i];
						}
						
						
						tn.addCost(this->cf.getCost(edge, edge2)
							* factor);	
		//				std::cout << "Edges 1cost+=" << cf.getCost(edge, edge2)* factor << " g1 edge:" << edge.from << "-" << edge.to << ":" << edge.elabel << " g2 edge:" << edge2.from << "-" << edge2.to << ":" << edge2.elabel << std::endl;
					}
					else {
						tn.addCost(this->cf.getEdgeCosts() * factor);  //	std::cout << "Edges 2cost+=" << cf.getEdgeCosts() * factor << std::endl;
					}
				}
				else { // deletion end2Index<0 means end2Index=-2
					tn.addCost(this->cf.getEdgeCosts() * factor);// std::cout << "Edges 3cost+=" << cf.getEdgeCosts() * factor << std::endl;
				}
			}
		}
		for (unsigned int j = 0; j < originalGraph1.size(); j++){
			if (unconnectedNode[j] == 0){
				if (tn.matching[j] >= 0){
					if (tn.originalGraph2.hasEdge(endIndex, tn.matching[j])) {	// if there is an edge in graph2, we deal with the edge substitution
						tn.addCost(this->cf.getEdgeCosts() * factor); 
			//			std::cout << " processEdges unconnected t-t': " << endIndex << "-" << tn.matching[j] << std::endl;
			//			std::cout << "Edges 2cost+=" << cf.getEdgeCosts() * factor << std::endl;
					}
				}
			}
		}

		// DUPLICATED CODE REFACTOR
	}
public:
	TreeNode() {

	}
	void show(){ 
		if ((unusedNodes1.size() == 0) || (unusedNodes2.size() == 0)) { std::cout << "At least one g no unmatched node" << " totalcost:" << (*this).estimatedTotalcost_unitCost() << " empty?u1,u2: " << this->unusedNodes1.empty() << this->unusedNodes2.empty() << " depth" << depth << std::endl; return; }
		std::cout << "First unused node in g1:" << unusedNodes1[0].v_id << " in g2:" << unusedNodes2[0].v_id 
		<< " empty?u1,u2: " << this->unusedNodes1.empty()<< this->unusedNodes2.empty()<<" cost"<<cost<<" depth"<<depth<<" totalcost:"<<(*this).estimatedTotalcost_unitCost() << std::endl;
	}

	void show_matching(){
		std::cout << "Matching:		";
		for (unsigned int i = 0; i < matching.size(); i++){
			std::cout << matching[i] << " ";
		}
		std::cout << std::endl<< "inverseMatching:	";
		for (unsigned int i = 0; i < inverseMatching.size(); i++){
			std::cout  << inverseMatching[i] << " ";
		}
		std::cout << std::endl;
	}

	bool operator < (const TreeNode& str) const
	{
		if (depth == str.depth)
			return (cost > str.cost);
		else
			return (depth < str.depth);
	}
	TreeNode(const PSS::Graph &g1,const PSS::Graph &g2, CostFunction cf2, double Factor) {
		unusedNodes1 = g1;
		unusedNodes1.Assign_vertexid();
		unusedNodes2 = g2;
		unusedNodes2.Assign_vertexid();
		originalGraph1 = g1;
		originalGraph1.Assign_vertexid();
		originalGraph2 = g2;
		originalGraph2.Assign_vertexid();
		cost = 0;
		cf = cf2;
		matching.resize(g1.size());
		inverseMatching.resize(g2.size());
		for (int i = 0; i < matching.size(); i++) {
			matching[i] = -1;
		}
		for (int i = 0; i < inverseMatching.size(); i++) {
			inverseMatching[i] = -1;
		}
		factor = Factor;
	}

	std::vector<TreeNode> generateSuccessors(double bound) {
		std::vector<TreeNode> successors;
		
		// all nodes of g2 are processed, the remaining nodes of g1 are deleted
		if (this->unusedNodes2.empty()) {
			TreeNode tn=(*this);
			int n = tn.unusedNodes1.size();
			int e = 0;
			for (unsigned int i = 0; i < n; i++){
				e += this->getNumberOfAdjacentEdges(tn.matching, originalGraph1, unusedNodes1[i].v_id);
				tn.matching[unusedNodes1[i].v_id] = -2; // -2 = deletion
			}
			tn.addCost(n * this->cf.getNodeCosts());
			tn.addCost(e * this->cf.getEdgeCosts() * factor);
			tn.unusedNodes1.clear();
	//		std::cout << "x u2 d depth " << tn.depth +1 << std::endl;// std::cin.get();
			if (this->cost <= bound){
				tn.depth++;
				successors.push_back(tn);
			}
		}
		else { // there are still nodes in g2 but no nodes in g1, the nodes of
			// g2 are inserted
			
			if (this->unusedNodes1.empty()) {
				
				TreeNode tn = (*this);
				int n = tn.unusedNodes2.size();
				int e = 0;
	//			std::cout << "x x d depth " << tn.depth + 1<< std::endl;
				for (unsigned int i = 0; i < n; i++){
					e += this->getNumberOfAdjacentEdges(tn.inverseMatching, originalGraph2, unusedNodes2[i].v_id); // Check here
					tn.inverseMatching[unusedNodes2[i].v_id] = -2; // -2 = deletion
				}
				tn.addCost(n * this->cf.getNodeCosts());
				tn.addCost(e * this->cf.getEdgeCosts() * factor);
				tn.unusedNodes2.clear();
				if (this->cost <= bound){
					tn.depth++;
					successors.push_back(tn);
				}
			}
			else { // there are nodes in both g1 and g2
				
				for (int i = 0; i < this->unusedNodes2.size(); i++) {	//In this loop, we deal with matching the first unused node in g1 to EACH of the unused node in g2
					
					TreeNode tn = (*this); 
					//Retrieves and removes the head(first element) of this list.
	//				PSS::Vertex start = tn.unusedNodes1.remove();
	//				PSS::Vertex end = tn.unusedNodes2.remove(i);
					tn.addCost(this->cf.getCost(tn.unusedNodes1[0], tn.unusedNodes2[i]));
	//				std::cout << "u1 u2 s " << unusedNodes1.Subgraph_Vertex_mapping.size()<< std::endl;
					tn.matching[unusedNodes1[0].v_id] = unusedNodes2[i].v_id;
					tn.inverseMatching[unusedNodes2[i].v_id] = unusedNodes1[0].v_id;
					
					this->processEdges(tn, unusedNodes1[0], unusedNodes1[0].v_id, unusedNodes2[i], unusedNodes2[i].v_id);
		//			std::cout << "u1 u2 s depth "<<tn.depth+1<<" " << tn.originalGraph1[tn.unusedNodes1.Subgraph_Vertex_mapping[0]].v_id << "-" <<
		//				tn.originalGraph2[tn.unusedNodes2.Subgraph_Vertex_mapping[i]].v_id <<" cost:"<<tn.cost<< std::endl;
					tn.unusedNodes1 = tn.unusedNodes1.remove_vertex(0);
					tn.unusedNodes2 = tn.unusedNodes2.remove_vertex(i);
					if (this->cost <= bound){
						tn.depth++;
						successors.push_back(tn); 
					}
				}			

				// deletion of a node from g_1 is also a valid successor
				TreeNode tn = (*this);	
				tn.matching[unusedNodes1[0].v_id] = -2;
				tn.addCost(this->cf.getNodeCosts());
				// find number of edges adjacent to node i
				int e = this->getNumberOfAdjacentEdges(tn.matching, originalGraph1, unusedNodes1[0].v_id);
				tn.addCost(this->cf.getEdgeCosts() * e
					* factor);
	//			std::cout << "u1 u2 d depth " << tn.depth + 1 << " " << tn.originalGraph1[tn.unusedNodes1.Subgraph_Vertex_mapping[0]].v_id << " cost:" << tn.cost <<  std::endl;
				tn.unusedNodes1 = tn.unusedNodes1.remove_vertex(0);
	//			Node deleted = tn.unusedNodes1.remove();
				if (this->cost <= bound){
					tn.depth++;
					successors.push_back(tn); 
				}
			}
		}
		return successors;
	}

	bool allNodesUsed() {
		if (unusedNodes1.empty() && unusedNodes2.empty()) {
			return true;
		}
		return false;
	}

	double getCost() {
		return cost;
	}

	double estimatedTotalcost_unitCost() {		
	// estimated total cost = exsiting cost of the matched part + lower bound of the cost for matching the unmatched nodes, 
	// in the case of unit cost, it is just the numbers of vertex and edge relabeling between the remaining part.
		if ((unusedNodes1.size() == 0) || (unusedNodes1.size() == 0)) return cost;
		
		std::vector<int> nodelable1, nodelable2, edgelable1, edgelable2;
		
		for (unsigned int i = 0; i < unusedNodes1.size(); i++){
			nodelable1.push_back(unusedNodes1[i].label);
			for (unsigned int j = 0; j < unusedNodes1[i].edge.size(); j++){
				edgelable1.push_back(unusedNodes1[i].edge[j].elabel);
			}
		}

		for (unsigned int i = 0; i < unusedNodes2.size(); i++){
			nodelable2.push_back(unusedNodes2[i].label);
			for (unsigned int j = 0; j < unusedNodes2[i].edge.size(); j++){
				edgelable2.push_back(unusedNodes2[i].edge[j].elabel);
			}
		}
		
		std::sort(nodelable1.begin(), nodelable1.end());
		std::sort(nodelable2.begin(), nodelable2.end());
		std::sort(edgelable1.begin(), edgelable1.end());
		std::sort(edgelable2.begin(), edgelable2.end());

		std::vector<int> v; v.resize(nodelable1.size()+nodelable2.size());
		std::vector<int> e; e.resize(edgelable1.size() + edgelable2.size());
		std::vector<int>::iterator it; 
		it = std::set_intersection(nodelable1.begin(), nodelable1.end(), nodelable2.begin(), nodelable2.end(), v.begin());
		int n1, n2;
		v.resize(it - v.begin()); 
		n1 = nodelable1.size() - v.size();
		n2 = nodelable2.size() - v.size();
		if (n1 < n2) n1 = n2; 
		it = std::set_intersection(edgelable1.begin(), edgelable1.end(), edgelable2.begin(), edgelable2.end(), e.begin());
		e.resize(it - e.begin());
		int e1, e2;
		e1 = edgelable1.size() - e.size();
		e2 = edgelable2.size() - e.size();
		if (e1 < e2) e1 = e2;


		return cost+n1+e1/2;
	}

	double estimatedTotalcost_arbitraryCost() {
		class MatrixGenerator cm;
		double estimatedFutureCost;
		Matrix<double> costMatrix = cm.getMatrix(unusedNodes1, unusedNodes2);		
		Hungarian _h(costMatrix);
		_h.max_to_min();
		estimatedFutureCost= _h.solve();
		return cost + estimatedFutureCost;
	}

};




}

#endif