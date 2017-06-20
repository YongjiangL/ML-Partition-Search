
namespace ull_ge{

class ulliso{
private:
	std::vector<int> m1;
	std::vector<int> m2;
	std::vector<std::vector<unsigned int>>	M;
	unsigned int d=0;
	bool flag = false;
	PSS::Graph S;
	PSS::Graph L;
public:
	ulliso(PSS::Graph &a, PSS::Graph &b){
		S = a;
		L = b;
		m1.resize(S.vertex_size());
		for (unsigned int i = 0; i < m1.size(); i++){ m1[i] = -1; }
		m2.resize(L.vertex_size());		
		M.resize(S.vertex_size());
		for (unsigned int i = 0; i < M.size(); i++){
			M.at(i).resize(L.vertex_size());
			for (unsigned int j = 0; j < M.at(i).size(); j++){
				if ((S[i].label == L[j].label) && (S[i].edge.size() <= L[j].edge.size())){
					M[i][j] = 1;
				}
				else
					M[i][j] = 0;
			}
		}
	}

	bool match(){				
		if (d == m1.size())	return true;

		for (unsigned int i = 0; i < m2.size(); i++){
			if ((m2[i] == 0)&&(refine(i)&&M[d][i])){
				d++;
				m2[i] = 1;
				m1[d-1] = i;					
				if (match()) return true;
				m1[d] = -1;
				d--;			
				m2[i] = 0;
			}
		}
		return false;
	}

	bool refine(int k){
		for (unsigned int j = 0; j < S[d].edge.size(); j++){
			int t = m1.at(S[d].edge[j].to);
			if (t != -1){
				if (!L.hasEdgeElabel(k, t, S[d].edge[j].elabel)) return false;
			}
		}
		return true;
	}

	void showmatch(){
		for (unsigned int i = 0; i < m1.size(); i++) std::cout << "U" << i << " match to V" << m1[i] << std::endl;
	}

};

class halfedge_ulliso{
private:
	std::vector<int> m1;
	std::vector<int> m2;
	std::vector<std::vector<unsigned int>>	M;
	unsigned int d = 0;
	bool flag = false;
	PSS::Graph S;
	PSS::Graph L;
	bool dummy=0;
public:
	halfedge_ulliso(PSS::Graph &a, PSS::Graph &b){
		S = a;
		L = b;
		if (S[S.size() - 1].label == -2) dummy = 1;	// Check if there is a dummy node in the graph
		if (dummy)
		m1.resize(S.vertex_size()-1);		// For half edge sub iso verification, we don't match dummy node. so here we have -1. Assuming S is half edge graph
		else
		m1.resize(S.vertex_size());

		for (unsigned int i = 0; i < m1.size(); i++){ m1[i] = -1; }
		m2.resize(L.vertex_size());
		
		if (dummy)
			M.resize(S.vertex_size() - 1);		// For half edge sub iso verification
		else
			M.resize(S.vertex_size());


		for (unsigned int i = 0; i < M.size(); i++){
			M.at(i).resize(L.vertex_size());
			for (unsigned int j = 0; j < M.at(i).size(); j++){
				if ((S[i].label == L[j].label) && (S[i].edge.size() <= L[j].edge.size())){
					M[i][j] = 1;
				}
				else
					M[i][j] = 0;
			}
		}
	}

	bool match(){//for (unsigned int i = 0; i < m1.size(); i++) std::cout << i << "-" << m1[i] << " "; std::cout << std::endl;
//		std::cout << "Depth:" << d <<" m1:"<<m1.size()<< std::endl;
		if (d == m1.size())	return true;
		
		for (unsigned int i = 0; i < m2.size(); i++){
			if ((m2[i] == 0) && (refine(i) && M[d][i])){
				d++;
				m2[i] = 1;
				m1[d - 1] = i;
				if (match()) return true;
				m1[d] = -1;
				d--;
				m2[i] = 0;
			}
		}
		return false;
	}

	bool refine(int k){
		for (unsigned int j = 0; j < S[d].edge.size(); j++){
			if (dummy){
				if (S[d].edge[j].to == S.size() - 1) continue;		// For dummy node

				int t = m1.at(S[d].edge[j].to);
				if (t != -1){
					if (!L.hasEdgeElabel(k, t, S[d].edge[j].elabel)) return false;
				}
			}
			else{
				int t = m1.at(S[d].edge[j].to);
				if (t != -1){
					if (!L.hasEdgeElabel(k, t, S[d].edge[j].elabel)) return false;
				}
			}
		}
		return true;
	}

	void showmatch(){
		for (unsigned int i = 0; i < m1.size(); i++) std::cout << "U" << i << " match to V" << m1[i] << std::endl;
	}

};




}