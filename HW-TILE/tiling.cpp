#ifndef BIPARTGRAPH_H
#define BIPARTGRAPH_H

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
	while (P[P.size() - 1] != s)
		P.push_back(prev[P[P.size() - 1]]);

	// Reverse shortest path
	for (int i = 0; i < P.size() / 2; ++i)
		swap(P[i], P[P.size() - 1 - i]);

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
		for (int i = 0; i < P.size() - 1; ++i)
		{
			--((*(resV.find(P[i])))->weights[P[i + 1]]);
			++((*(resV.find(P[i + 1])))->weights[P[i]]);
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


class BiPartGraph
{
private:

	//Stores the number of vertices, used to determine the first condition in has_tiling
	int numVertices;
	//Stores the set with all the vertices
	unordered_set<Vertex*> totalCheckers;
	//Stores the set of black checkers in the bipartite graph
	unordered_set<Vertex*> blackCheckers;
	//Stores the set of red checkers in the bipartite graph
	unordered_set<Vertex*> redCheckers;
	//Stores all the vertices in the graph with their coordinates
	unordered_map<Vertex*, pair<int, int>> vertexDictionary;

	//The vertices will be used for max flow and be computed for perfect matching
	Vertex *source;
	Vertex *sink;



	void addNeighbors()
	{
		string dummy = "";
		for (auto i : blackCheckers)
		{
			pair<int, int> up, down, left, right;
			//Set up coordinates
			up.first = vertexDictionary.at(i).first - 1;
			up.second = vertexDictionary.at(i).second;
			//Set down coordinates
			down.first = vertexDictionary.at(i).first + 1;
			down.second = vertexDictionary.at(i).second;
			//Set left coordinates
			left.first = vertexDictionary.at(i).first;
			left.second = vertexDictionary.at(i).second - 1;
			//Set right coordinates
			right.first = vertexDictionary.at(i).first;
			right.second = vertexDictionary.at(i).second + 1;

			//Step 1: Try to search for a neighbor up
			for (auto k : redCheckers)
			{
				if (vertexDictionary.at(k) == up)
				{
					i->neighs.insert(k);
					i->weights[k] = 1;
				}
				else if (vertexDictionary.at(k) == down)
				{
					i->neighs.insert(k);
					i->weights[k] = 1;
				}
				else if (vertexDictionary.at(k) == left)
				{
					i->neighs.insert(k);
					i->weights[k] = 1;
				}
				else if (vertexDictionary.at(k) == right)
				{
					i->neighs.insert(k);
					i->weights[k] = 1;
				}
			}
		}
	}

	//This function connects the source to all the blac checkers
	void setSource()
	{
		for (auto i : blackCheckers)
		{
			source->neighs.insert(i);
			//Sets the max flow to 1
			source->weights[i] = 1;
		}
		totalCheckers.insert(source);
	}

	//This function connects all the red checkers to the sink
	void setSink()
	{
		for (auto i : redCheckers)
		{
			i->neighs.insert(sink);
			//Sets the max flow to 1
			i->weights[sink] = 1;
		}
		totalCheckers.insert(sink);
	}


public:

	//Default constructor
	BiPartGraph()
	{
		source = new Vertex();
		sink = new Vertex();
		numVertices = 0;
	}

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
				blackCheckers.insert(baby);
				totalCheckers.insert(baby);
				vertexDictionary[baby] = bPair;
				column++;
			}
			else
			{
				Vertex * baby = new Vertex();
				pair<int, int> bPair;
				bPair.first = row;
				bPair.second = column;
				redCheckers.insert(baby);
				totalCheckers.insert(baby);
				vertexDictionary[baby] = bPair;
				column++;
			}
		}
		addNeighbors();
		setSource();
		setSink();
	}
	Vertex* GetSource()
	{
		return source;
	}
	Vertex* GetSink()
	{
		return sink;
	}

	int getFlow()
	{
		return max_flow(source, sink, totalCheckers);
	}

	int getB()
	{
		int counter = 0;
		for (auto i : blackCheckers)
		{
			counter++;
		}
		return counter;
	}

	///Helper method to display variables
	void displayFlow()
	{
		int counter = 1;

		cout << "Source" << endl;
		cout << "Neighbors: ";
		for (auto x : source->neighs)
			cout << vertexDictionary[x].first << "," << vertexDictionary[x].second << " :: ";
		cout << endl;

		for (auto i : vertexDictionary)
		{
			cout << counter << ": " << i.second.first << "," << i.second.second << endl;
			cout << "Neighbors: ";
			for (auto k : i.first->neighs)
			{
				cout << vertexDictionary[k].first << "," << vertexDictionary[k].second << " :: ";
			}
			counter++;
			cout << endl;
		}

		cout << "Sink" << endl;
		cout << "Neighbors: ";
		for (auto y : sink->neighs)
			cout << vertexDictionary[y].first << "," << vertexDictionary[y].second << " :: ";
		cout << endl;

		system("pause");
	}
};

bool has_tiling(string floor)
{
	string modFloor = "";
	int row = 0;
	int column = 0;
	int flow, numB;
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

	BiPartGraph CheckerBoard;

	CheckerBoard.constructGraph(modFloor);

	if (CheckerBoard.isValid() == false)
	{
		return false;
	}

	Vertex *s = CheckerBoard.GetSource();
	Vertex *t = CheckerBoard.GetSink();

	//max flow
	flow = CheckerBoard.getFlow();
	numB = CheckerBoard.getB();

	//augmented path
	//cout << modFloor << endl;
	//CheckerBoard.displayFlow();
	
	if (flow == numB)
		return true;
	else
		return false;
}


#endif // !BIPARTGRAPH_H