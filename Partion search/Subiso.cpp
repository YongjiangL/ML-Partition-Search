#include "stdafx.h"
#include "Subiso.h"
#include "Database.h"

namespace ull{

UllState::UllState(PSS::Graph *g1, PSS::Graph *g2)
{
	n1 = g1->size();
	n2 = g2->size();

	core_len = 0;

	core_1 = new node_id[n1];
	core_2 = new node_id[n2];
	M = new byte *[n1];
	if (!core_1 || !core_2 || !M)
	{
		using namespace std;
		cout << "Out of memory" << endl;
		exit(-1);
	}

	int i, j;

	for (i = 0; i<n1; i++)
	{
		M[i] = new byte[n2];
		if (!M[i])
		{
			using namespace std;
			cout << "Out of memory" << endl;
			exit(-1);
		}
	}
	
	for (i = 0; i<n1; i++)
	{
		core_1[i] = NULL_NODE;
	}
	for (i = 0; i<n2; i++)
	{
		core_2[i] = NULL_NODE;
	}
	for (i = 0; i<n1; i++)
	for (j = 0; j < n2; j++)
		M[i][j] = (g1->at(i).edge.size() == g2->at(i).edge.size() &&
//		g1->OutEdgeCount(i) == g2->OutEdgeCount(j)) &&
//g1->CompatibleNode(g1->GetNodeAttr(i), g2->GetNodeAttr(j)) 
		g1->at(i).label == g2->at(i).label )
		?	1 : 0;

}

void UllState::GetCoreSet(node_id c1[], node_id c2[])
{
	int i, j;
	for (i = 0, j = 0; i<n1; i++)
	if (core_1[i] != NULL_NODE)
	{
		c1[j] = i;
		c2[j] = core_1[i];
		j++;
	}
}

bool UllState::NextPair(node_id *pn1, node_id *pn2,
	node_id prev_n1, node_id prev_n2)
{
	if (prev_n1 == NULL_NODE)
	{
		prev_n1 = core_len;
		prev_n2 = 0;
	}
	else if (prev_n2 == NULL_NODE)
		prev_n2 = 0;
	else
		prev_n2++;

	if (prev_n2 >= n2)
	{
		prev_n1++;
		prev_n2 = 0;
	}

	if (prev_n1 != core_len)
		return false;
	while (prev_n2<n2 && M[prev_n1][prev_n2] == 0)
		prev_n2++;
	if (prev_n2<n2)
	{
		*pn1 = prev_n1;
		*pn2 = prev_n2;
		return true;
	}
	else
		return false;
}

bool UllState::IsFeasiblePair(node_id node1, node_id node2)
{
	assert(node1<n1);
	assert(node2<n2);

	return M[node1][node2] != 0;
}

UllState* UllState::Clone()
{
	return new UllState(*this);
}

void UllState::AddPair(node_id node1, node_id node2)
{
	assert(node1<n1);
	assert(node2<n2);
	assert(core_len<n1);
	assert(core_len<n2);

	core_1[node1] = node2;
	core_2[node2] = node1;

	core_len++;

	int k;

	for (k = core_len; k<n1; k++)
		M[k][node2] = 0;

	refine();
}

UllState::~UllState()
{
	delete[] core_1;
	delete[] core_2;
	int i;
	for (i = 0; i<n1; i++)
	if (M[i])
		delete[] M[i];
	delete[] M;
}

}