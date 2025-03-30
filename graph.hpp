// namespace graph(non-directed graph) with adjacency list (arrays)

#include <iostream>
#include <stdexcept>
#include <new>

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
        Edge() : src(0), dest(0), weight(0) {}
        Edge(int src, int dest, int weight) : src(src), dest(dest), weight(weight) {}
        ~Edge() {}
    };

    // Graph class
    class Graph
    {
    private:
        // function to allocate memory for the adjacency list
        Node *getAdjacencyList(int val, int weight, Node *head)
        {
            // allocate memory for the new node
            Node *newNode = new (std::nothrow) Node;
            if (!newNode)
                throw std::bad_alloc();

            newNode->val = val;
            newNode->cost = weight;
            // point the new node to the head of the adjacency list
            newNode->next = head;
            return newNode;
        }
        // constant number of vertices
        const int numVertices;
        // adjacency list
        Node **head;

    public:
        // constructors
        Graph(int numVertices) : numVertices(numVertices)
        {
            // check if numVertices is valid
            if (numVertices <= 0)
                throw std::invalid_argument("Number of vertices must be greater than 0!");

            // allocate memory for the adjacency list
            this->head = new (std::nothrow) Node *[numVertices]();
            if (!this->head)
                throw std::bad_alloc();

            for (int i = 0; i < numVertices; i++)
            {
                head[i] = nullptr;
            }
        }

        Graph(Edge edges[], int n, int numVertices);

        // destructor
        ~Graph();

        // check if edge exists
        bool hasEdge(int src, int dest) const;

        // add an edge to the graph
        void addEdge(int src, int dest, int weight = 1);

        // remove an edge from the graph
        void removeEdge(int src, int dest);

        int getNumVertices() const;

        // get the adjacency list
        Node *getHead(int src) const;

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