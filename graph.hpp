// namespace graph(non-directed graph) with adjacency list (arrays)

namespace graph
{

    struct Node
    {
        int val, cost;
        Node *next;
    };

    struct Edge
    {
        int src, dest, weight;
    };

    // Graph class
    class Graph
    {
        // function to allocate memory for the adjacency list
        Node *getAdjacencyList(int val, int weight, Node *head)
        {
            Node *newNode = new Node;
            newNode->val = val;
            newNode->cost = weight;
            // point the new node to the head of the adjacency list
            newNode->next = head;
            return newNode;
        }

    private:
        // number of vertices
        int numVertices;

    public:
        // adjacency list
        Node **head;

        // constructor
        Graph(Edge edges[], int n, int numVertices);

        // destructor
        ~Graph();

        // check if edge exists
        bool hasEdge(int src, int dest) const;

        // add an edge to the graph
        void addEdge(int src, int dest, int weight = 1);

        // remove an edge from the graph
        void removeEdge(int src, int dest);

        // print the graph
        void print_graph(Node *ptr, int i);
    };

    class Algorithms
    {
    public:
        // return a the sub tree that the root is start
        static Graph bfs(const Graph &graph, int start);

        // return a tree or forest that the root is start its
        // will contain all the vertices that are reachable from the start and tree edges
        static Graph dfs(const Graph &graph, int start);

        // return a tree or forest that the root is start and
        // the weight of the edge is the shortest path
        static Graph dijkstra(const Graph &graph, int start);

        // return a spanning tree of the graph
        static Graph prim(const Graph &graph);

        // return a spanning tree of the graph
        static Graph kruskal(const Graph &graph);
    };
}