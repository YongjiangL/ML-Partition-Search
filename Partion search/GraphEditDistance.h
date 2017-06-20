#ifndef GraphEditDistance_h
#define GraphEditDistance_h
#include <iostream>
#include <vector>
#include "Database.h"
#include "TreeNode.h"
#include <map>
#include <string>
#include <vector>
#include <limits>
#include "Hungarian.h"
#include <algorithm>    // std::sort
#include <queue> 
namespace GED{



class EditDistance {
private:
	bool undirected;	//whether or not the edges are directed
public: 
	EditDistance(int undirected) {
		this->undirected = undirected;
	}
	EditDistance() {
		this->undirected = 1;
	}
	/**
	*
	* @return the exact edit distance between graph @param g1
	* and graph @param g2 using the cost function @param cf
	* s is the maximum number of open paths used in beam-search
	*/
	double getEditDistance(PSS::Graph g1, PSS::Graph g2, CostFunction cf, double bound,int c) {
		std::priority_queue<TreeNode> open;
	// if the edges are undirected
	// all of the edge operations have to be multiplied by 0.5
	// since all edge operations are performed twice 
		double factor = 1.0;
		if (this->undirected){
		factor = 1;
		}
		TreeNode minGed;
		TreeNode start(g1, g2, cf, factor); 
		open.push(start);
		std::vector<TreeNode> successors;
		double MinimumEditDistance = std::numeric_limits<double>::infinity();
		
		while (!open.empty()){
			
			TreeNode u = open.top();
//			std::cout << std::endl << "Tree loop. open.size=" << open.size() << std::endl; u.show(); u.show_matching(); //std::cin.get();
			open.pop(); 
			if (u.estimatedTotalcost_unitCost()>MinimumEditDistance){		// Prune unpromising branch
				continue;
			}
			if (c){			
				if (u.estimatedTotalcost_arbitraryCost() >MinimumEditDistance){		// Prune unpromising branch
				continue;}
			}

		
			if (u.allNodesUsed()){
				if (u.getCost() < MinimumEditDistance)
				{
					MinimumEditDistance = u.getCost();
					minGed = u;
				}
				if (MinimumEditDistance < bound)
					bound = MinimumEditDistance;
				continue;
			}
			// generates all successors of node u in the search tree
			// and add them to open
			
			successors = u.generateSuccessors(bound);
//			std::cout << "successors.size=" << successors.size() << std::endl;
			for (unsigned int i = 0; i < successors.size(); i++){
				open.push(successors[i]);
			}
			successors.clear();
//			while (open.size() > s){
	//			open.pop_back();
		//	}
		}
//		minGed.show_matching();
		return MinimumEditDistance;
	}

	/**
	*
	* @return the approximated edit distance between graph @param g1
	* and graph @param g2 according to the @param matching using the cost function @param cf
	*/
	double getEditDistance(PSS::Graph g1, PSS::Graph g2, Matrix<int> matching, CostFunction cf) {
		double factor = 1.0;
		if (undirected){
			factor = 0.5;
		}
		Matrix<int> edgesOfG1 = g1.ToAdjacencyMatrix();
		Matrix<int> edgesOfG2 = g2.ToAdjacencyMatrix();
		vector<double> impliedEdgeCosts(matching.columns(), 0);
		double ed = 0;
		int maxCostMatch = -1;
		double maxEdgeCost = 0;
		for (int i = 0; i < matching.columns(); i++){
			if (matching.m_matrix[i][0] < g1.size()){
				if (matching.m_matrix[i][1] < g2.size()){
					// i-th node substitution with node from g2 with index matching[i][1]
					ed += cf.getCost(g1[matching.m_matrix[i][0]], g2[matching.m_matrix[i][1]]);
					double impliedEdgeCosts = 0;
					
				}
			}
		}
	}
};


class BipartiteMatching {
private:
	double assignmentCost = 0;
public:
	BipartiteMatching() {}
	double getAssignmentCost() {
		return assignmentCost;
	}
	/**
	* @return the optimal matching according to the @param costMatrix
	* the matching actually used is defined in the string "matching"
	*/
	Matrix<int> getMatching(Matrix<double> costMatrix) {
		Matrix<int>	assignment;
//		HungarianAlgorithm ha = new HungarianAlgorithm();
//		assignment = ha.hgAlgorithm(costMatrix);
		Hungarian _h(costMatrix);
		_h.max_to_min();
//		assignment = _h.solve();
//		cout << _h.solve() << endl;
		//			this.assignmentCost = ha.hgAlgorithmOnlyCost(costMatrix);
		//			RiesenMunkres ha = new RiesenMunkres(costMatrix);
		//			assignment = ha.getMinCostMatching();
		return assignment;
	}
};

class GraphMatching {
private:
	/**
	* the source and target graph actually to be matched (temp is for temporarily swappings)
	*/
	PSS::Graph sourceGraph, targetGraph, temp;
	double distance;
	/**
	* whether the edges of the graphs are undirected (=1) or directed (=0)
	*/
	int undirected;
	/**
	* the cost function to be applied
	*/
	CostFunction costFunction;
	/**
	* computes an optimal bipartite matching of local graph structures
	*/
	BipartiteMatching bipartiteMatching;
	/**
	* computes the approximated or exact graph edit distance
	*/
	EditDistance editDistance;
	/**
	* the matching procedure defined via GUI or properties file
	* possible choices are 'Hungarian', 'VJ' (VolgenantJonker)
	* 'AStar' (exact tree search) or 'Beam' (approximation based on tree-search)
	*/
	std::string matching;
	/**
	* generates the cost matrix whereon the optimal bipartite matching can
	* be computed
	*/
	MatrixGenerator matrixGenerator;
	std::string adj; // best or worst or none
public:
	/**
	* the matching procedure
	*/
	GraphMatching(PSS::Graph sourceGraph, PSS::Graph targetGraph) {

	// the cost matrix used for bipartite matchings
	Matrix<double> costMatrix;
	// distance value d
	double d = -1;
	// if both graphs are empty the distance is zero and no computations have to be carried out!
	if (sourceGraph.size()<1 && targetGraph.size()<1){
		d = 0;
	}
	else{
		// calculate the approximated or exact edit-distance using tree search algorithms 
		// AStar: number of open paths during search is unlimited (s=infty)
		// Beam: number of open paths during search is limited to s
			if (matching == "AStar"){
				d = editDistance.getEditDistance(
				sourceGraph, targetGraph, costFunction, std::numeric_limits<double>::infinity(),0);
			}
			else {
			// approximation of graph edit distances via bipartite matching
			// in order to get determinant edit costs between two graphs
				if (sourceGraph.size()<targetGraph.size()){
					temp = sourceGraph;
					sourceGraph = targetGraph;
					targetGraph = temp;
				}
				// compute the matching using Hungarian or VolgenantJonker (defined in String matching)
		//		int[][] matching = null;
				// generate the cost-matrix between the local substructures of the source and target graphs
				costMatrix = matrixGenerator.getMatrix(sourceGraph, targetGraph);
				Matrix<int> bi_matching = bipartiteMatching.getMatching(costMatrix);
				d = editDistance.getEditDistance(sourceGraph,targetGraph, bi_matching, costFunction);
			}
	}
	}


};

}

#endif