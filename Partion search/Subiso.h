#ifndef Subiso_h
#define Subiso_h
#ifndef BYTE_DEFINED
#define BYTE_DEFINED
typedef unsigned char byte;
#endif
#include "Database.h"

namespace ull{

	typedef unsigned short node_id;
	const node_id NULL_NODE = 0xFFFF;

	/*
	class State
	{

	public:
		virtual ~State() {}
		virtual PSS::Graph *GetGraph1() = 0;
		virtual PSS::Graph *GetGraph2() = 0;
		virtual bool NextPair(node_id *pn1, node_id *pn2,
			node_id prev_n1 = NULL_NODE, node_id prev_n2 = NULL_NODE) = 0;
		virtual bool IsFeasiblePair(node_id n1, node_id n2) = 0;
		virtual void AddPair(node_id n1, node_id n2) = 0;
		virtual bool IsGoal() = 0;
		virtual bool IsDead() = 0;
		virtual int CoreLen() = 0;
		virtual void GetCoreSet(node_id c1[], node_id c2[]) = 0;
		virtual State *Clone() = 0;  // Changed clone to Clone for uniformity

		virtual void BackTrack() { };
	};
*/
	class UllState 
	{
	private:
		int core_len;
		node_id *core_1;
		node_id *core_2;
		PSS::Graph *g1, *g2;
		int n1, n2;
		byte **M;   // Matrix encoding the compatibility of the nodes

		void refine(){}

	public:
		UllState(PSS::Graph *g1, PSS::Graph *g2);
		bool IsGoal() { return core_len == n1 && core_len == n2; };
		int CoreLen() { return core_len; }
		void GetCoreSet(node_id c1[], node_id c2[]);		
		bool IsDead() {
			if (n1 != n2) return true;
			for (int i = core_len; i<n1; i++)
			{
				for (int j = 0; j<n2; j++)
				if (M[i][j] != 0) goto next_row;
				return true;
			next_row:;
			}
			return false;
		};		
		bool NextPair(node_id *pn1, node_id *pn2,
			node_id prev_n1 = NULL_NODE, node_id prev_n2 = NULL_NODE); 
		bool IsFeasiblePair(node_id n1, node_id n2);
		UllState *Clone();
		void AddPair(node_id n1, node_id n2);
		~UllState();
		void BackTrack() { };
		/*
		UllState(const UllState &Ullstate);
		
		PSS::Graph *GetGraph1() { return g1; }
		PSS::Graph *GetGraph2() { return g2; }

		
		
		*/

		
		
		
	};
	
	static bool match(int *pn, node_id c1[], node_id c2[], UllState *s)
	{bool found = false;
		if (s->IsGoal())
		{
			*pn = s->CoreLen();
			s->GetCoreSet(c1, c2);
			return true;
		}

		if (s->IsDead())
			return false;

		node_id n1 = NULL_NODE, n2 = NULL_NODE;
		
		while (!found && s->NextPair(&n1, &n2, n1, n2))
		{
			if (s->IsFeasiblePair(n1, n2))
			{
				UllState *s1 = s->Clone();
				s1->AddPair(n1, n2);
				found = match(pn, c1, c2, s1);
//				s1->BackTrack();
				delete s1;
			}
		}
		return found;
	}


}

#endif