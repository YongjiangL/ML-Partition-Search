#ifndef _HUNGARIAN_H_
#define _HUNGARIAN_H_

#include "Matrix.h"

#define INF 100000000        // just infinity

using namespace std;

class Hungarian		// input: a matrix output: the MAXIMUM assignment. If we need minimum assignment, run max_to_min as preprosessing
{
public:
	Hungarian(void);
	Hungarian(Matrix<double> &);
	~Hungarian(void);
	
public:
	void set(Matrix<double> &);
	double solve(void);								//O(n^3) 
	double solve(Matrix<double> &);
	double incre_solve(Matrix<double> &); // given new parts.      Warning, finding minimum assignment is not supported in incremental method now, need slight modification
	Matrix<double> cost; // cost matrix	// when set to finding minimum assignment, Mij will be changed to max-Mij
	Matrix<double> original_cost; // original cost matrix
	void max_to_min(void);				// Change the setting to find minimum assignment by set Mij to max-Mij
private:
	void init_labels(void);
	void init_labels(int, int);
	void update_labels(void);
	void add_to_tree(int, int);
	void augment(void);

private:
	void initAll(void);
	void freeAll(void);

private:
	
	int rows, columns;   // the row and column number before resize
	int n, max_match;    // n left nodes and n right nodes
	double *lx, *ly;     // labels of X and Y parts
	int *xy;             // xy[x] - vertex that is matched with x,
	int *yx;             // yx[y] - vertex that is matched with y
	bool *S, *T;         // sets S and T in algorithm
	double *slack;       // as in the algorithm description
	int *slackx;         // slackx[y] such a vertex, that
						 // l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
	int *prev;           // array for memorizing alternating paths
};

#endif
