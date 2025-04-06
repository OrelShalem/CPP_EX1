// orel8155@gmail.com
/**
 * @file graph.hpp
 * @brief Header file for the Graph implementation and graph algorithms
 *
 * This file defines the Graph class and Algorithms class for graph operations.
 * It includes data structures for graph representation (nodes and edges) and
 * declarations for various graph algorithms such as BFS, DFS, Dijkstra, Prim,
 * and Kruskal.
 */

#pragma once
#include <iostream>
#include <stdexcept>
#include <new>

namespace graph
{
    /**
     * @brief Node structure for adjacency list representation
     *
     * Represents a node in the adjacency list of a graph.
     * Each node contains the vertex value, the edge weight (cost),
     * and a pointer to the next node in the adjacency list.
     */
    struct Node
    {
        int val;    ///< The vertex that this node connects to
        int cost;   ///< The weight/cost of the edge
        Node *next; ///< Pointer to the next adjacent node
    };

    /**
     * @brief Edge structure for representing graph edges
     *
     * Represents an edge in the graph with source vertex,
     * destination vertex, and weight of the edge.
     */
    struct Edge
    {
        int src;    ///< Source vertex of the edge
        int dest;   ///< Destination vertex of the edge
        int weight; ///< Weight/cost of the edge
    };

    /**
     * @brief Graph class for representing and manipulating undirected graphs
     *
     * This class implements an undirected weighted graph using adjacency lists.
     * It provides methods for adding and removing edges, checking if edges exist,
     * and accessing the adjacency list of any vertex.
     */
    class Graph
    {
    private:
        /**
         * @brief Allocates memory for a new node in the adjacency list.
         *
         * @param val The target vertex for this adjacency list entry
         * @param weight The weight of the edge
         * @param head The current head of the adjacency list
         * @return Node* Pointer to the newly created node
         * @throws std::bad_alloc if memory allocation fails
         */
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
        Node **head;     ///< Array of adjacency lists, one for each vertex
        int numVertices; ///< Number of vertices in the graph
        int numEdges;    ///< Number of edges in the graph

    public:
        /**
         * @brief Constructs a graph with a specified number of vertices
         *
         * Creates an empty graph with the given number of vertices and no edges.
         *
         * @param numVertices Number of vertices in the graph
         * @throws std::invalid_argument If number of vertices is non-positive
         */
        explicit Graph(int numVertices) : numVertices(numVertices), numEdges(0)
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

        /**
         * @brief Constructs a graph with predefined edges
         *
         * Creates a graph from an array of edges, adding each edge to the graph.
         *
         * @param edges Array of edges to add to the graph
         * @param n Number of edges in the array
         * @param numVertices Number of vertices in the graph
         * @throws std::invalid_argument If vertex numbers are invalid
         */
        Graph(Edge edges[], int n, int numVertices);

        /**
         * @brief Destructor for the Graph class
         *
         * Frees all dynamically allocated memory for adjacency lists.
         */
        ~Graph();

        /**
         * @brief Checks if an edge exists between two vertices
         *
         * @param src Source vertex
         * @param dest Destination vertex
         * @return true If the edge exists
         * @return false If the edge doesn't exist
         * @throws std::invalid_argument If vertex numbers are invalid
         */
        bool hasEdge(int src, int dest) const;

        /**
         * @brief Adds an edge between two vertices
         *
         * Adds an undirected edge between the source and destination vertices
         * with the specified weight.
         *
         * @param src Source vertex
         * @param dest Destination vertex
         * @param weight Weight of the edge (default: 1)
         * @throws std::invalid_argument If vertex numbers are invalid, if it's a self-loop,
         *         or if the edge already exists
         */
        void addEdge(int src, int dest, int weight = 1);

        /**
         * @brief Removes an edge between two vertices
         *
         * Removes the undirected edge between the source and destination vertices
         * if it exists.
         *
         * @param src Source vertex
         * @param dest Destination vertex
         * @throws std::invalid_argument If vertex numbers are invalid or if the edge doesn't exist
         */
        void removeEdge(int src, int dest);

        /**
         * @brief Gets the number of vertices in the graph
         *
         * @return int Number of vertices
         */
        int getNumVertices() const;

        /**
         * @brief Gets the number of edges in the graph
         *
         * @return int Number of edges
         */
        int getNumEdges() const;

        /**
         * @brief Gets the head of the adjacency list for a vertex
         *
         * @param src Vertex number
         * @return Node* Pointer to the head of the adjacency list
         * @throws std::invalid_argument If the vertex number is invalid
         */
        Node *getHead(int src) const;

        /**
         * @brief Prints the adjacency list for a vertex
         *
         * @param ptr Head of the adjacency list to print
         * @param i Vertex number
         * @throws std::invalid_argument If the vertex number is invalid
         */
        void print_graph(Node *ptr, int i);
    };

    /**
     * @brief Creates a new adjacency list node
     *
     * Helper function to create a new node for the adjacency list.
     *
     * @param val Vertex value
     * @param weight Edge weight
     * @param next Pointer to the next node
     * @return Node* New node for the adjacency list
     */
    Node *getAdjacencyList(int val, int weight, Node *next);

    /**
     * @brief Class containing static methods for graph algorithms
     *
     * This class provides implementations of various graph algorithms
     * such as BFS, DFS, Dijkstra, Prim, and Kruskal.
     */
    class Algorithms
    {
    public:
        /**
         * @brief Performs a Breadth-First Search on the graph
         *
         * Constructs and returns a tree representing the BFS traversal
         * starting from the specified vertex.
         *
         * @param graph Input graph
         * @param start Starting vertex
         * @return Graph A tree representing the BFS traversal
         * @throws std::invalid_argument If the start vertex is invalid
         */
        static Graph bfs(const Graph &graph, int start);

        /**
         * @brief Performs a Depth-First Search on the graph
         *
         * Constructs and returns a tree representing the DFS traversal
         * starting from the specified vertex.
         *
         * @param graph Input graph
         * @param start Starting vertex
         * @return Graph A tree representing the DFS traversal
         * @throws std::invalid_argument If the start vertex is invalid
         */
        static Graph dfs(const Graph &graph, int start);

        /**
         * @brief Performs Dijkstra's algorithm for finding shortest paths
         *
         * Constructs and returns a tree representing the shortest paths
         * from the start vertex to all other vertices.
         *
         * @param graph Input graph
         * @param start Starting vertex
         * @return Graph A tree representing the shortest paths
         * @throws std::invalid_argument If the start vertex is invalid
         */
        static Graph dijkstra(const Graph &graph, int start);

        /**
         * @brief Performs Prim's algorithm for finding a Minimum Spanning Tree
         *
         * Constructs and returns a tree representing the Minimum Spanning Tree (MST)
         * of the input graph.
         *
         * @param graph Input graph
         * @return Graph A tree representing the MST
         * @throws std::runtime_error If the graph doesn't have enough edges (n-1)
         */
        static Graph prim(const Graph &graph);

        /**
         * @brief Performs Kruskal's algorithm for finding a Minimum Spanning Tree
         *
         * Constructs and returns a tree representing the Minimum Spanning Tree (MST)
         * of the input graph.
         *
         * @param graph Input graph
         * @return Graph A tree representing the MST
         * @throws std::runtime_error If the graph doesn't have enough edges (n-1)
         */
        static Graph kruskal(const Graph &graph);

        // /**
        //  * @brief Checks if the graph is connected
        //  *
        //  * Determines whether there is a path between any two vertices in the graph.
        //  *
        //  * @param graph Input graph
        //  * @return true If the graph is connected
        //  * @return false If the graph has disconnected components
        //  */
        // static bool isConnected(const Graph &graph);
    };
}