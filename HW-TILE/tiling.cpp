

#include "tiling.h"
#include "vertex.h"

using namespace std;


// Finds a (shortest according to edge length) augmenting path
// from s to t in a graph with vertex set V.
// Returns whether there is an augmenting path.
bool augmenting_path(Vertex* s, Vertex* t, unordered_set<Vertex*> V, vector<Vertex*> &P)
{
        // Check that s and t aren't nullptr
        if (s == nullptr || t == nullptr)
	{
		cerr << "augmenting_path() was passed nullptr s or t." << endl;
		abort();
	}

        // Check that s and t are in the graph
        if (V.find(s) == V.end() || V.find(t) == V.end())
	{
		cerr << "augmenting_path() was passed s or t not in V." << endl;
		abort();
	}

	// Check that every vertex has valid neighs/weights.
	for (Vertex* v : V)
		for (Vertex* vn : v->neighs)
			if (v->weights.find(vn) == v->weights.end())
			{
				cerr << "augmenting_path() was passed invalid vertex." << endl;
				abort();
			}

        // Since augmenting paths should have the fewest edges,
	// not the minimum weight, run BFS.
	queue<Vertex*> Q;
	Q.push(s);

	unordered_set<Vertex*> R;
	R.clear(); 
	R.insert(s);

	unordered_map<Vertex*, Vertex*> prev;

	while (!Q.empty())
	{
		Vertex* cur = Q.front();
		Q.pop();

		for (Vertex* nei : cur->neighs)
		{
			// Must have positive edge weight
			if (cur->weights[nei] == 0)
				continue;

			if (R.find(nei) == R.end())
			{
				Q.push(nei);
				R.insert(nei);
				prev[nei] = cur; 
			}
		}
	}      

        // If BFS never reached t
        if (R.find(t) == R.end())
                return false;

        // Reconstruct shortest path backwards
        P.clear();
        P.push_back(t);
        while (P[P.size()-1] != s)
                P.push_back(prev[P[P.size()-1]]);

        // Reverse shortest path
        for (int i = 0; i < P.size()/2; ++i)
		swap(P[i], P[P.size()-1-i]);

        return true;
}

// Returns the maximum flow from s to t in a weighted graph with vertex set V.
// Assumes all edge weights are non-negative.
int max_flow(Vertex* s, Vertex* t, unordered_set<Vertex*> V)
{
	// If s or t is invalid.
        if (s == nullptr || t == nullptr)
	{
		cerr << "max_flow() was passed nullptr s or t." << endl;
		abort(); 
	}

	// If s or t is not in the vertex set.
        if (V.find(s) == V.end() || V.find(t) == V.end())
	{
		cerr << "max_flow() was passed s or t not in V." << endl;
		abort(); 
	}

	// Check that every vertex has valid neighs/weights.
	for (Vertex* v : V)
		for (Vertex* vn : v->neighs)
			if (v->weights.find(vn) == v->weights.end())
			{
				cerr << "max_flow() was passed invalid vertex." << endl;
				abort();
			}

        // Create a deep copy of V to use as the residual graph
        unordered_set<Vertex*> resV;
        unordered_map<Vertex*, Vertex*> C; // Maps vertices in V to copies in resV
        for (Vertex* vp : V)
        {
                Vertex* rp = new Vertex;
                resV.insert(rp);
                C[vp] = rp;
        }
        for (Vertex* vp : V)
                for (Vertex* np : vp->neighs)
                {
                        C[vp]->neighs.insert(C[np]);
                        C[vp]->weights[C[np]] = vp->weights[np];
                }
	// Add any missing necessary "back" edges. 
        for (Vertex* vp : V)
                for (Vertex* np : vp->neighs)
		{
			if (C[np]->neighs.find(C[vp]) == C[np]->neighs.end())
			{
				C[np]->neighs.insert(C[vp]);
				C[np]->weights[C[vp]] = 0;
			}
		}

        // Run Edmonds-Karp
        while (true)
        {
                // Find an augmenting path
                vector<Vertex*> P;
                if (!augmenting_path(C[s], C[t], resV, P))
                        break;  
                // Update residual graph
                for (int i = 0; i < P.size()-1; ++i)
                {
                        --((*(resV.find(P[i])))->weights[P[i+1]]);
                        ++((*(resV.find(P[i+1])))->weights[P[i]]);
                }
        }

        // Compute actual flow amount
        int flow = 0;
        for (Vertex* snp : C[s]->neighs)
                flow += 1 - C[s]->weights[snp];

        // Delete residual graph
        for (Vertex* vp : resV)
                delete vp;

        return flow;
}

inline size_t key(int i, int j) { return (size_t)i << 32 | (unsigned int)j; }

class BiPartGraph
{
private:

	//Stores the number of vertices, used to determine the first condition in has_tiling
	int numVertices;
	//Stores the set of black checkers in the bipartite graph
	unordered_map<size_t, Vertex*> blackCheckers;
	//Stores the set of red checkers in the bipartite graph
	unordered_map<size_t, Vertex*> redCheckers;
	//Stores all the vertices in the graph
	unordered_map<Vertex*, pair<int, int>> vertexDictionary;

	void addNeighbors()
	{
		string dummy = "";
		for (auto i : blackCheckers)
		{
			pair<int, int> up, down, left, right;
			//Set up coordinates
			up.first = vertexDictionary.at(i.second).first - 1;
			up.second = vertexDictionary.at(i.second).second;
			//Set down coordinates
			down.first = vertexDictionary.at(i.second).first + 1;
			down.second = vertexDictionary.at(i.second).second;
			//Set left coordinates
			left.first = vertexDictionary.at(i.second).first;
			left.second = vertexDictionary.at(i.second).second - 1;
			//Set right coordinates
			right.first = vertexDictionary.at(i.second).first;
			right.second = vertexDictionary.at(i.second).second + 1;

			//Step 1: Try to search for a neighbor up
			try
			{
				Vertex* neighbor = redCheckers.at(key(up.first, up.second));
				i.second->neighs.insert(neighbor);
				i.second->weights[neighbor] = 1;
			}
			catch (out_of_range oor)
			{
				dummy = "";
			}
			//Step 2: Try to search for a neighbor down
			try
			{
				Vertex* neighbor = redCheckers.at(key(down.first, down.second));
				i.second->neighs.insert(neighbor);
				i.second->weights[neighbor] = 1;
			}
			catch(out_of_range oor1)
			{
				dummy = "";
			}
			//Step 3: Try to search for a neighbor to the left
			try
			{
				Vertex* neighbor = redCheckers.at(key(left.first, left.second));
				i.second->neighs.insert(neighbor);
				i.second->weights[neighbor] = 1;
			}
			catch (out_of_range oor2)
			{
				dummy = "";
			}
			//Step 4: Try to search for a neighbor to the right
			try
			{
				Vertex* neighbor = redCheckers.at(key(right.first, right.second));
				i.second->neighs.insert(neighbor);
				i.second->weights[neighbor] = 1;
			}
			catch (out_of_range oor3)
			{
				dummy = "";
			}
		}
	}

public:

	//This function returns false if the two sets do not have the same number of elements
	bool isValid()
	{
		if (blackCheckers.size() == redCheckers.size())
			return true;
		else
			return false;
	}

	void constructGraph(string floor)
	{
		int row = 0;
		int column = 0;

		for (int i = 0; i < floor.length(); i++)
		{
			if (floor[i] == '#')
				column++;
			else if (floor[i] == '\n')
			{
				row++;
				column = 0;
			}
			else if (floor[i] == 'b')
			{
				Vertex * baby = new Vertex();
				pair<int, int> bPair;
				bPair.first = row;
				bPair.second = column;
				blackCheckers[key(row, column)] = baby;
				vertexDictionary[baby] = bPair;
				column++;
			}
			else
			{
				Vertex * baby = new Vertex();
				pair<int, int> bPair;
				bPair.first = row;
				bPair.second = column;
				redCheckers[key(row, column)] = baby;
				vertexDictionary[baby] = bPair;
				column++;
			}
		}
		//addNeighbors();
	}
};

bool has_tiling(string floor)
{
	string modFloor = "";
	int row = 0;
	int column = 0;
	bool firstElmnt = false;
	bool startsEven;

	//This loops transform the string into another string formated as a checkers board
	for (int i = 0; i < floor.length(); i++)
	{
		if (floor[i] == '#')
		{
			modFloor += floor[i];
			column++;
		}
		else if (floor[i] == '\n')
		{
			row++;
			column = 0;
			modFloor += floor[i];
		}
		else if (floor[i] == ' ')
		{
			if (!firstElmnt)
			{
				if (column % 2 == 0)
					startsEven = true;
				else
					startsEven = false;
				modFloor += 'b';
				firstElmnt = true;
			}
			else if (startsEven)
			{
				if (row % 2 != 0)
				{
					if (column % 2 == 0)
						modFloor += 'b';
					else
						modFloor += 'r';
				}
				else
				{
					if (column % 2 != 0)
						modFloor += 'b';
					else
						modFloor += 'r';
				}
			}
			else
			{
				if (row % 2 != 0)
				{
					if (column % 2 == 0)
						modFloor += 'r';
					else
						modFloor += 'b';
				}
				else
				{
					if (column % 2 != 0)
						modFloor += 'r';
					else
						modFloor += 'b';
				}
			}
			column++;
		}
	}


	//cout << modFloor << endl;
	//system("pause");
    return false;
}




