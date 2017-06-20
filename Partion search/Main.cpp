// Main.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Database.h"
#include "iostream"
#include "Ullman.h"
#include "GraphEditDistance.h"
#include "timer.h"
using namespace std;

int k = 2;		// number of partitions for each graph in the database
int threshold = 2;   //The GED threshold
string graphDataBaseFileName = "dataset/data.graphs.2";
string queryGraphFileName = "dataset/query.graphs.1";
string outputFileName = "./output/results";
int profileFilter = 1;

void test_Ullman();
void test_Canonicalcode();
void test_GED();


int _tmain(int argc, _TCHAR* argv[])
{	
	cout << "data graphs file: " << graphDataBaseFileName << endl;
	cout << "query graphs file: " << queryGraphFileName << endl;
	cout << "output file: " << outputFileName << endl;

	PSS::Database D;
	D.IndexingMethods_push(1);
	D.IndexingMethods_push(3);
	D.Init(graphDataBaseFileName, queryGraphFileName, outputFileName, k, threshold);
	D.Build_index(graphDataBaseFileName, k, outputFileName);
	D.Query_processing(k, threshold, outputFileName, profileFilter);

	cout << "Done!" << endl;
	return 0;
}

void test_Ullman() {
	PSS::Database D;
	D.Init(graphDataBaseFileName, queryGraphFileName, outputFileName, k, threshold);
	ull_ge::halfedge_ulliso is(D.Graphs[0], D.Graphs[1]);
	cout << "Matched?" << is.match() << endl;
	is.showmatch();
}

void test_Canonicalcode() {
	PSS::DFSCode minCode, minCode2;
	class TSINGHUA_CLIPSE_UTIL::TimeRecorder time;
	int time1 = 0, time2 = 1;
	float t_time;
	PSS::Database D;
	D.IndexingMethods_push(1);
	D.Init(graphDataBaseFileName, queryGraphFileName, outputFileName, k, threshold);
	cout << "Number of graphs:" << D.Graphs.size() << endl;
	std::cout << "---------------------------------------------" << endl;
	for (int i = 0; i <D.Graphs.size() && (i+1) < D.Graphs.size(); i = i + 2){
		time.check(); time1++; time2++;
		minCode = D.CanonicalCode(D.Graphs[i]);
		std::cout << "GraphID: " << D.Graphs[i].ID << endl;
		time.check();
		t_time = time.diffTime(time1++, time2++);
		std::cout << "Time for building code " << i << "  size:" << D.Graphs[i].size() << " is:" << t_time << endl << minCode.toString() << endl;

		time.check(); time1++; time2++;
		minCode2 = D.CanonicalCode(D.Graphs[i + 1]);

		time.check();
		t_time = time.diffTime(time1++, time2++);
		std::cout << "Time for building code " << i << "  size:" << D.Graphs[i + 1].size() << " is:" << t_time << endl << minCode2.toString() << endl;
		std::cout << "---------------------------------------------" << endl;
	}
}

void test_GED() {
	cout << "testing GED of two graphs" << endl;
	PSS::Database D;
	D.Init(graphDataBaseFileName, queryGraphFileName, outputFileName, k, threshold);

	PSS::Graph a = D.Graphs[1];
	PSS::Graph b = D.Graphs[2];
	a.printGraph("./output/a"); 
	cout << "------------ " << endl;
	b.printGraph("./output/b");

	class TSINGHUA_CLIPSE_UTIL::TimeRecorder time;
	int time1 = 0, time2 = 1;
	float t_time;
	time.check(); time1++; time2++;

	GED::CostFunction cf;
	GED::EditDistance ed;
	cout << "GED:" << ed.getEditDistance(a, b, cf, 4, 0) << endl;
	time.check();
	t_time = time.diffTime(time1++, time2++);
	std::cout << "GED test time:" << t_time << endl;
}
