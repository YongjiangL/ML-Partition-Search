#include "stdafx.h"
#include "Hungarian.h"
#include <iostream>
#include <fstream>
#include <string>

#define max(a,b) (a>b)?a:b
#define min(a,b) (a<b)?a:b

using namespace std;

Hungarian::Hungarian(void)
{
	initAll();
}

Hungarian::Hungarian(Matrix<double> &cost_matrix)
{
	initAll();
	
	this->rows = cost_matrix.rows();
	this->columns = cost_matrix.columns();
	
	this->cost = cost_matrix;
	this->n = cost_matrix.maxsize();
	this->cost.resize(this->n, this->n);
	this->original_cost = this->cost;
}


void Hungarian::max_to_min(void) {
	double max = 0;
	for (int i = 0; i < cost.rows(); i++)
	for (int j = 0; j < cost.columns(); j++)
	if (cost.m_matrix[i][j]>max)
		max = cost.m_matrix[i][j];
	for (int i = 0; i < cost.rows(); i++)
	for (int j = 0; j < cost.columns(); j++)
		cost.m_matrix[i][j] = max - cost.m_matrix[i][j];	
}

Hungarian::~Hungarian(void)
{
	freeAll();
}

void Hungarian::set(Matrix<double> &cost_matrix)
{
	this->rows = cost_matrix.rows();
	this->columns = cost_matrix.columns();

	this->cost = cost_matrix;
	this->n = cost_matrix.maxsize();
	this->cost.resize(this->n, this->n);
	this->original_cost = this->cost;
}

void Hungarian::initAll(void)
{
	rows = 0;
	columns = 0;
	n = 0;
	max_match = 0;
	lx = NULL;
	ly = NULL;
	xy = NULL;
	yx = NULL;
	S = NULL;
	T = NULL;
	prev = NULL;
	slack = NULL;
	slackx = NULL;
}

void Hungarian::freeAll(void)
{
	if (lx != NULL)
		delete[] lx;
	if (ly != NULL)
		delete[] ly;
	if (xy != NULL)
		delete[] xy;
	if (yx != NULL)
		delete[] yx;
	if (S != NULL)
		delete[] S;
	if (T != NULL)
		delete[] T;
	if (prev != NULL)
		delete[] prev;
	if (slack != NULL)
		delete[] slack;
	if (slackx != NULL)
		delete[] slackx;

	lx = NULL;
	ly = NULL;
	xy = NULL;
	yx = NULL;
	S = NULL;
	T = NULL;
	prev = NULL;
	slack = NULL;
	slackx = NULL;
}

void Hungarian::init_labels(void)
{
	lx = new double[this->n];
	ly = new double[this->n];
	for (int i = 0; i < this->n; i++)
	{
		lx[i] = 0.0;
		ly[i] = 0.0;
	}
	for (int x = 0; x < this->n; x++)
	for (int y = 0; y < this->n; y++)
	{
		lx[x] = max(lx[x], cost(x, y));
//		cout << "cost " <<x<<"-"<<y<<" is "<<cost(x,y) << endl;
	}
			
}

void Hungarian::init_labels(int old_rows, int old_cols)
{
    double *nlx = this->lx;
	double *nly = this->ly;
	this->lx = new double[n];
	this->ly = new double[n];

	int *nxy = this->xy;
	int *nyx = this->yx;
	this->xy = new int[n];
	this->yx = new int[n];

	for (int i = 0, j = 0 ; i < n && j < n; i++, j++)
	{
		if (i < old_rows)
		{
			this->lx[i] = nlx[i];
			this->xy[i] = nxy[i];
		}
		else
		{
			this->lx[i] = 0.0;
			this->xy[i] = -1;			
		}

		if (j < old_cols)
		{
			this->ly[j] = nly[j];
			this->yx[j] = nyx[j];
		}
		else
		{
			this->ly[j] = 0.0;
			this->yx[j] = -1;
		}
	}
	for (int y = old_cols; y < n; y++)
	{
		for (int x = 0; x < old_rows; x++)
			this->ly[y] = max(ly[y], cost(x, y) - lx[x]);

		this->ly[y] = max(ly[y], cost(y, y));
	}
	for (int x = old_rows; x < n; x++)
	{
		for (int y = 0; y < n; y++)
			this->lx[x] = max(lx[x], cost(x, y) - ly[y]); 
	}

	delete[] nlx;
	nlx = NULL;
	delete[] nly;
	nly = NULL;
	delete[] nxy;
	nxy = NULL;
	delete[] nyx;
	nyx = NULL;
}

void Hungarian::update_labels(void)
{
    int x, y;
	double delta = INF;                //init delta as infinity
    for (y = 0; y < n; y++)            //calculate delta using slack
        if (!T[y])
            delta = min(delta, slack[y]);
    for (x = 0; x < n; x++)            //update X labels
        if (S[x]) lx[x] -= delta;
    for (y = 0; y < n; y++)            //update Y labels
        if (T[y]) ly[y] += delta;
    for (y = 0; y < n; y++)            //update slack array
        if (!T[y])
            slack[y] -= delta;
}

void Hungarian::add_to_tree(int x, int prevx) 
//x - current vertex,prevx - vertex from X before x in the alternating path,
//so we add edges (prevx, xy[x]), (xy[x], x)
{
    S[x] = true;                    //add x to S
    prev[x] = prevx;                //we need this when augmenting
    for (int y = 0; y < n; y++)     //update slacks, because we add new vertex to S
        if (lx[x] + ly[y] - cost(x, y) < slack[y])
        {
            slack[y] = lx[x] + ly[y] - cost(x, y);
            slackx[y] = x;
        }
}

void Hungarian::augment(void)          //main function of the algorithm
{
    if (max_match == n) return;        //check wether matching is already perfect
    int x, y, root;                    //just counters and root vertex
    int *q = new int[n];
	int wr = 0, rd = 0;          //q - queue for bfs, wr,rd - write and read
                                       //pos in queue
	S = new bool[n];
	T = new bool[n];
	prev = new int[n];
	for (int i = 0; i < n; i++)
	{
		S[i] = false;
		T[i] = false;
		prev[i] = -1;
	}
    for (x = 0; x < n; x++)            //finding root of the tree
        if (xy[x] == -1)
        {
			q[wr++] = root = x;//	cout << endl<<"wr is "<<wr-1<<" q[wr-1]"<<q[wr-1] << endl;
            prev[x] = -2;
			S[x] = true;		//	cout  << "x " << x << " S[x] " << S[x] << endl;
            break;
        }
		
    for (y = 0; y < n; y++)            //initializing slack array
    {
        slack[y] = lx[root] + ly[y] - cost(root, y);
		slackx[y] = root;		//cout << endl << "root "<<root<<" lx[root] " << lx[root] << " ly[y] " << ly[y] << " cost(root, y) " << cost(root, y) << endl;
    }


	while (true)                                                        //main cycle
    {
        while (rd < wr)                                                 //building tree with bfs cycle
        {
            x = q[rd++];                                                //current vertex from X part
            for (y = 0; y < n; y++)                                     //iterate through all edges in equality graph
                if (cost(x, y) == lx[x] + ly[y] &&  !T[y])
                {
                    if (yx[y] == -1) break;                             //an exposed vertex in Y found, so
                                                                        //augmenting path exists!
                    T[y] = true;                                        //else just add y to T,
                    q[wr++] = yx[y];                                    //add vertex yx[y], which is matched
                                                                        //with y, to the queue
                    add_to_tree(yx[y], x);                              //add edges (x,y) and (y,yx[y]) to the tree
                }
            if (y < n) break;                                           //augmenting path found!
        }
        if (y < n) break;                                               //augmenting path found!

        update_labels();                                                //augmenting path not found, so improve labeling
        wr = rd = 0;                
        for (y = 0; y < n; y++)        
        //in this cycle we add edges that were added to the equality graph as a
        //result of improving the labeling, we add edge (slackx[y], y) to the tree if
        //and only if !T[y] &&  slack[y] == 0, also with this edge we add another one
        //(y, yx[y]) or augment the matching, if y was exposed
            if (!T[y] &&  slack[y] == 0)
            {
                if (yx[y] == -1)                                        //exposed vertex in Y found - augmenting path exists!
                {
                    x = slackx[y];
                    break;
                }
                else
                {
                    T[y] = true;                                        //else just add y to T,
                    if (!S[yx[y]])    
                    {
                        q[wr++] = yx[y];                                //add vertex yx[y], which is matched with
                                                                        //y, to the queue
                        add_to_tree(yx[y], slackx[y]);                  //and add edges (x,y) and (y,
                                                                        //yx[y]) to the tree
                    }
                }
            }
        if (y < n) break;                                               //augmenting path found!
    }

    if (y < n)                                                          //we found augmenting path!
    {
        max_match++;                                                    //increment matching
        //in this cycle we inverse edges along augmenting path
        for (int cx = x, cy = y, ty; cx != -2; cx = prev[cx], cy = ty)
        {
            ty = xy[cx];
            yx[cy] = cx;
            xy[cx] = cy;
        }
		delete[] S;
		delete[] T;
		delete[] prev;
		S = NULL;
		T = NULL;
		prev = NULL;
        augment();                          //recall function, go to step 1 of the algorithm
    }

	delete [] q;
}//end of augment() function

double Hungarian::solve(Matrix<double> &cost_matrix)
{
	this->rows = cost_matrix.rows();
	this->columns = cost_matrix.columns();

	this->cost = cost_matrix;
	this->n = cost_matrix.maxsize();
	this->cost.resize(this->n, this->n);
	this->original_cost = this->cost;
	double ret = 0;
    max_match = 0;                       //number of vertices in current matching
	xy = new int[n];	
    yx = new int[n];
	for (int i = 0 ; i < n; i++)
	{
		xy[i] = -1;
		yx[i] = -1;
	}
	slack = new double[n];
	slackx = new int[n];
	for (int i = 0; i < n; i++)
	{
		slack[i] = 0.0;
		slackx[i] = -1;
	}
    init_labels();                    //step 0
    augment();                        //steps 1-3
	for (int row = 0; row < original_cost.rows(); row++){
		for (int col = 0; col < original_cost.columns(); col++){
			if (xy[row] == col){
				ret += original_cost(row, col);
				cost_matrix(row,col) = 0;
			}
			else
				cost_matrix(row,col) = -1;
		}
	}

	return ret;
}

double Hungarian::solve(void)
{
	double ret = 0;
    max_match = 0;                       //number of vertices in current matching
	xy = new int[n];	
    yx = new int[n];
	for (int i = 0 ; i < n; i++)
	{
		xy[i] = -1;
		yx[i] = -1;
	}
	slack = new double[n];
	slackx = new int[n];
	for (int i = 0; i < n; i++)
	{
		slack[i] = 0.0;
		slackx[i] = -1;
	}
    init_labels();                    //step 0
    augment();                        //steps 1-3
    for (int x = 0; x < n; x++)       //forming answer there
	{
//		cout << "x=" << x << " xy[x]=" << xy[x] << " cost=" << original_cost(x, xy[x]) << endl;
        ret += original_cost(x, xy[x]);
		if (original_cost(x, xy[x]) == 0)
		{
			yx[xy[x]] = -1;
			xy[x] = -1;
			--max_match;
		}
	}

	return ret;
}

double Hungarian::incre_solve(Matrix<double> &cost_matrix)
{
	if (max_match == 0)
	{
		cout << "wrong call the function for incremental pass!" << endl;
		exit(-1);
	}

	int old_rows = this->rows;
	int old_cols = this->columns;
	this->rows = cost_matrix.rows();
	this->columns = cost_matrix.columns();

	this->cost = cost_matrix;
	this->n = cost_matrix.maxsize();
	this->cost.resize(this->n, this->n);

	double ret = 0;                      //weight of the optimal matching
	init_labels(old_rows, old_cols);

	if (slack != NULL)
	{
		delete[] slack;
		slack = NULL;
	}
	if (slackx != NULL)
	{
		delete[] slackx;
		slackx = NULL;
	}
	slack = new double[n];
	slackx = new int[n];
	for (int i = 0; i < n; i++)
	{
		slack[i] = 0.0;
		slackx[i] = -1;
	}
	augment();                           //steps 1-3

	for ( int row = 0; row < cost.rows(); row++ ){
			for ( int col = 0 ; col < cost.columns() ; col++ ){
				if (xy[row] == col){
					ret += cost(row,col);
					cost_matrix(row,col) = 0;
				}
				else
					cost_matrix(row,col) = -1;
			}
	}

    return ret;
}
