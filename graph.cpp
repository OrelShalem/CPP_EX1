// orel8155@gmail.com
#include "graph.hpp"
#include <iostream>
#include <climits>
#include "Queue.hpp"
#include "PriorityQueue.hpp"
#include "UnionFind.hpp"

namespace graph
{
    /**
     * @brief Constructor for creating a graph from an array of edges
     *
     * This constructor initializes a graph with the given number of vertices
     * and adds the edges provided in the array. It handles self-loops and
     * duplicate edges by issuing warnings and ignoring them.
     *
     * @param edges Array of edges to be added to the graph
     * @param n Number of edges in the array
     * @param numVertices Number of vertices in the graph
     */
    Graph::Graph(Edge edges[], int n, int numVertices) : numVertices(numVertices), numEdges(0)
    {
        this->head = new Node *[numVertices]();

        for (int i = 0; i < numVertices; i++)
        {
            head[i] = nullptr;
        }

        for (int i = 0; i < n; i++)
        {
            int src = edges[i].src;
            int dest = edges[i].dest;
            int weight = edges[i].weight;

            try
            {
                // Use addEdge instead of directly adding to the list
                addEdge(src, dest, weight);
            }
            catch (const std::invalid_argument &e)
            {
                // If addEdge throws an exception, print a warning and continue
                std::cout << "Warning: " << e.what() << std::endl;
            }
        }
    }

    /**
     * @brief Destructor for the Graph class
     *
     * Frees all dynamically allocated memory for the adjacency lists
     * and the head array itself.
     */
    Graph::~Graph()
    {
        for (int i = 0; i < numVertices; i++)
        {
            Node *current = head[i];
            while (current != nullptr)
            {
                Node *temp = current;
                current = current->next;
                delete temp;
            }
        }
        delete[] head;
    }

    /**
     * @brief Checks if an edge exists between two vertices
     *
     * Validates the vertices and then traverses the adjacency list
     * of the source vertex to find a node with the destination vertex value.
     *
     * @param src Source vertex
     * @param dest Destination vertex
     * @return true If the edge exists
     * @return false If the edge doesn't exist
     * @throws std::invalid_argument If either vertex is out of range
     */
    bool Graph::hasEdge(int src, int dest) const
    {
        // check if src and dest are valid vertices
        if (src >= getNumVertices() || dest >= getNumVertices() || src < 0 || dest < 0)
            // throw an error
            throw std::invalid_argument("Invalid vertex numbers!");

        // check if edge exists
        Node *current = head[src];
        while (current != nullptr)
        {
            if (current->val == dest)
                return true;
            current = current->next;
        }
        return false;
    }

    /**
     * @brief Adds an edge between two vertices
     *
     * Validates the vertices, checks for self-loops and duplicate edges,
     * then adds the edge to the adjacency lists of both vertices (since this is an undirected graph).
     * Increments the edge count.
     *
     * @param src Source vertex
     * @param dest Destination vertex
     * @param weight Weight of the edge (default: 1)
     * @throws std::invalid_argument If vertices are invalid, if it's a self-loop, or if the edge already exists
     */
    void Graph::addEdge(int src, int dest, int weight)
    {
        // check if src and dest are valid vertices
        if (src >= getNumVertices() || dest >= getNumVertices() || src < 0 || dest < 0)
        {
            throw std::invalid_argument("Invalid vertex numbers!");
        }

        // check for self-loop
        if (src == dest)
        {
            throw std::invalid_argument("Cannot add self-loop to vertex " + std::to_string(src));
        }

        // check for duplicate edge
        if (hasEdge(src, dest))
        {
            throw std::invalid_argument("Edge already exists between " + std::to_string(src) + " and " + std::to_string(dest));
        }

        Node *newNode = getAdjacencyList(dest, weight, head[src]);
        head[src] = newNode;

        // undirected graph
        newNode = getAdjacencyList(src, weight, head[dest]);
        head[dest] = newNode;

        numEdges++;
    }

    /**
     * @brief Removes an edge between two vertices
     *
     * Validates the vertices, checks if the edge exists, then
     * removes the edge from the adjacency lists of both vertices.
     * Decrements the edge count.
     *
     * @param src Source vertex
     * @param dest Destination vertex
     * @throws std::invalid_argument If vertices are invalid or if the edge doesn't exist
     */
    void Graph::removeEdge(int src, int dest)
    {
        // check if src and dest are valid vertices
        if (src >= getNumVertices() || dest >= getNumVertices() || src < 0 || dest < 0)
        {
            throw std::invalid_argument("Invalid vertex numbers!");
        }

        // check if edge exists
        if (!hasEdge(src, dest))
        {
            throw std::invalid_argument("Edge does not exist between " + std::to_string(src) + " and " + std::to_string(dest));
        }

        Node *current = head[src];
        Node *prev = nullptr;

        // Search for the edge to remove
        while (current != nullptr)
        {
            if (current->val == dest)
            {
                if (prev == nullptr)
                    head[src] = current->next;
                else
                    prev->next = current->next;
                delete current;
                break;
            }
            prev = current;
            current = current->next;
        }

        current = head[dest];
        prev = nullptr;
        while (current != nullptr)
        {
            if (current->val == src)
            {
                if (prev == nullptr)
                    head[dest] = current->next;
                else
                    prev->next = current->next;
                delete current;
                break;
            }
            prev = current;
            current = current->next;
        }

        numEdges--;
    }

    /**
     * @brief Gets the number of vertices in the graph
     *
     * @return int Number of vertices
     */
    int Graph::getNumVertices() const
    {
        return numVertices;
    }

    /**
     * @brief Gets the number of edges in the graph
     *
     * @return int Number of edges
     */
    int Graph::getNumEdges() const
    {
        return numEdges;
    }

    /**
     * @brief Gets the head of the adjacency list for a vertex
     *
     * @param src Source vertex
     * @return Node* Pointer to the head of the adjacency list
     * @throws std::invalid_argument If the vertex is invalid
     */
    Node *Graph::getHead(int src) const
    {
        if (src >= getNumVertices() || src < 0)
            throw std::invalid_argument("Invalid vertex number!");
        return head[src];
    }

    /**
     * @brief Prints the adjacency list for a vertex
     *
     * Displays the vertex and all its neighbors along with edge weights.
     *
     * @param ptr Pointer to the head of the adjacency list
     * @param i Vertex index
     * @throws std::invalid_argument If the vertex is invalid
     */
    void Graph::print_graph(Node *ptr, int i)
    {
        if (i >= getNumVertices() || i < 0)
            throw std::invalid_argument("Invalid vertex number!");

        // print the graph
        std::cout << "Vertex " << i << " -> ";
        if (ptr == nullptr)
        {
            std::cout << "No adjacent vertices";
        }
        else
        {
            std::cout << "Adjacent vertices: ";
            while (ptr != nullptr)
            {
                std::cout << "\n    To vertex " << ptr->val
                          << " (weight: " << ptr->cost << ")";
                ptr = ptr->next;
            }
        }
        std::cout << "\n"
                  << std::endl;
    }

    /**
     * @brief Breadth-First Search algorithm implementation
     *
     * Performs a BFS traversal starting from the given vertex and
     * constructs a tree representing the BFS traversal.
     *
     * @param graph The input graph
     * @param start The starting vertex
     * @return Graph A tree representing the BFS traversal
     * @throws std::invalid_argument If the start vertex is invalid
     */
    Graph Algorithms::bfs(const Graph &graph, int start)
    {
        int numVertices = graph.getNumVertices();

        // Validate the start vertex
        if (start < 0 || start >= numVertices)
        {
            throw std::invalid_argument("Invalid start vertex!");
        }

        // Create an empty result graph
        Graph result(numVertices);

        // Array to mark visited vertices
        bool *visited = new bool[numVertices]();

        // Queue to store discovered but not yet processed vertices
        Queue q(numVertices);

        // Mark the start vertex as visited
        visited[start] = true;

        // Add the start vertex to the queue
        q.enqueue(start);

        // While the queue is not empty
        while (!q.isEmpty())
        {
            // Remove a vertex from the queue
            int u = q.dequeue();

            // Get the adjacency list of the vertex
            Node *neighbor = graph.getHead(u);

            // For each neighbor
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // If the neighbor has not been visited
                if (!visited[v])
                {
                    // Mark the neighbor as visited
                    visited[v] = true;

                    // Add an edge to the result graph
                    result.addEdge(u, v, weight);

                    // Add the neighbor to the queue
                    q.enqueue(v);
                }

                // Move to the next neighbor
                neighbor = neighbor->next;
            }
        }

        // Free memory
        delete[] visited;

        return result;
    }

    /**
     * @brief Utility function for DFS traversal
     *
     * Helper function for the DFS algorithm that recursively visits
     * vertices and constructs the DFS tree.
     *
     * @param graph The input graph
     * @param u Current vertex being visited
     * @param visited Array to track visited vertices
     * @param result Graph representing the DFS tree
     * @param parent Parent vertex of the current vertex
     */
    void dfsUtil(const Graph &graph, int u, bool visited[], Graph &result, int parent)
    {
        // Mark the current vertex as visited
        visited[u] = true;

        // Get the adjacency list of the vertex
        Node *neighbor = graph.getHead(u);

        // For each neighbor
        while (neighbor != nullptr)
        {
            int v = neighbor->val;
            int weight = neighbor->cost;

            // If the neighbor has not been visited
            if (!visited[v])
            {
                // Add an edge to the result graph
                result.addEdge(u, v, weight);

                // Recursive call to the neighbor - passing current vertex as parent
                dfsUtil(graph, v, visited, result, u);
            }
            // Prevent adding back edge to parent (for undirected graphs)
            // If the neighbor is the parent of the current vertex, skip it
            else if (v != parent)
            {
                // Here we can add additional logic if we want to handle back edges
                // For example: cycle detection
            }

            // Move to the next neighbor
            neighbor = neighbor->next;
        }
    }

    /**
     * @brief Depth-First Search algorithm implementation
     *
     * Performs a DFS traversal starting from the given vertex and
     * constructs a tree representing the DFS traversal.
     *
     * @param graph The input graph
     * @param start The starting vertex
     * @return Graph A tree representing the DFS traversal
     * @throws std::invalid_argument If the start vertex is invalid
     */
    Graph Algorithms::dfs(const Graph &graph, int start)
    {
        int numVertices = graph.getNumVertices();

        // Validate the start vertex
        if (start < 0 || start >= numVertices)
        {
            throw std::invalid_argument("Invalid start vertex!");
        }

        // Create an empty result graph
        Graph result(numVertices);

        // Array to mark visited vertices
        bool *visited = new bool[numVertices]();

        // Perform DFS from the start vertex
        dfsUtil(graph, start, visited, result, -1);

        // Free memory
        delete[] visited;

        return result;
    }

    /**
     * @brief Dijkstra's algorithm for finding shortest paths
     *
     * Implements Dijkstra's algorithm to find the shortest paths from
     * a given source vertex to all other vertices in the graph.
     * Returns a tree representing these shortest paths.
     *
     * @param graph The input graph
     * @param start The starting vertex
     * @return Graph A tree representing the shortest paths
     * @throws std::invalid_argument If the start vertex is invalid
     */
    Graph Algorithms::dijkstra(const Graph &graph, int start)
    {
        int numVertices = graph.getNumVertices();

        // Validate the start vertex
        if (start < 0 || start >= numVertices)
        {
            throw std::invalid_argument("Invalid start vertex!");
        }

        // Create an empty result graph
        Graph result(numVertices);

        // Array of shortest distances from start vertex
        int *distance = new int[numVertices];

        // Array to mark vertices that are finalized
        bool *finalized = new bool[numVertices]();

        // Array to store the parent of each vertex in the shortest path
        int *parent = new int[numVertices];

        // Initialize arrays
        for (int i = 0; i < numVertices; i++)
        {
            distance[i] = INT_MAX; // Infinity
            parent[i] = -1;        // No parent
        }

        // Distance from start vertex to itself is 0
        distance[start] = 0;

        // Priority queue to find the next vertex with minimum distance
        PriorityQueue pq;

        // Add the start vertex to the priority queue
        pq.enqueue(start, 0);

        // While the queue is not empty
        while (!pq.isEmpty())
        {
            // Remove the vertex with minimum distance
            int u = pq.peek();
            pq.dequeue();

            // If we've already finalized this vertex, continue
            if (finalized[u])
                continue;

            // Mark this vertex as finalized
            finalized[u] = true;

            // If it has a parent, add the edge to the result graph
            if (parent[u] != -1)
            {
                // The weight is the distance difference, not the original edge weight!
                result.addEdge(parent[u], u, distance[u] - distance[parent[u]]);
            }

            // Process all neighbors of the current vertex
            Node *neighbor = graph.getHead(u);
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // If the neighbor is not finalized and there's a shorter path through u
                if (!finalized[v] && distance[u] != INT_MAX &&
                    distance[u] + weight < distance[v])
                {
                    // Update the distance
                    distance[v] = distance[u] + weight;

                    // Update the parent
                    parent[v] = u;

                    // Add to the priority queue
                    pq.enqueue(v, distance[v]);
                }

                // Move to the next neighbor
                neighbor = neighbor->next;
            }
        }

        // Free memory
        delete[] distance;
        delete[] finalized;
        delete[] parent;

        return result;
    }

    /**
     * @brief Prim's algorithm for finding Minimum Spanning Tree
     *
     * Implements Prim's algorithm to find a Minimum Spanning Tree (MST)
     * of the given graph. Requires at least n-1 edges to form an MST.
     *
     * @param graph The input graph
     * @return Graph A tree representing the MST
     * @throws std::runtime_error If the graph doesn't have at least n-1 edges
     */
    Graph Algorithms::prim(const Graph &graph)
    {
        // throw an error if the graph has not has at least n-1 edges
        if (graph.getNumEdges() < graph.getNumVertices() - 1)
        {
            throw std::runtime_error("Graph does not have at least n-1 edges!");
        }

        int numVertices = graph.getNumVertices();

        // Create an empty result graph
        Graph result(numVertices);

        // Array to mark vertices already included in MST
        bool *inMST = new bool[numVertices]();

        // Array of minimum weights to connect each vertex to the tree
        int *key = new int[numVertices];

        // Array of parents of each vertex in the spanning tree
        int *parent = new int[numVertices];

        // Initialize arrays
        for (int i = 0; i < numVertices; i++)
        {
            key[i] = INT_MAX; // Infinity
            parent[i] = -1;   // No parent
        }

        // Start with the first vertex (0)
        key[0] = 0; // Weight to the start vertex is 0

        // Priority queue to find the next vertex with minimum weight
        PriorityQueue pq;

        // Add the start vertex to the priority queue
        pq.enqueue(0, 0);

        // While we haven't added all vertices to the tree
        while (!pq.isEmpty())
        {
            // Remove the vertex with minimum weight
            int u = pq.peek();
            pq.dequeue();

            // If we've already added this vertex to the tree, continue
            if (inMST[u])
                continue;

            // Add the vertex to the tree
            inMST[u] = true;

            // If it has a parent, add the edge to the result graph
            if (parent[u] != -1)
            {
                result.addEdge(u, parent[u], key[u]);
            }

            // Process all neighbors of the current vertex
            Node *neighbor = graph.getHead(u);
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // If the neighbor is not in the tree and there's a lighter edge
                if (!inMST[v] && weight < key[v])
                {
                    // Update the weight
                    key[v] = weight;

                    // Update the parent
                    parent[v] = u;

                    // Add to the priority queue
                    pq.enqueue(v, key[v]);
                }

                // Move to the next neighbor
                neighbor = neighbor->next;
            }
        }

        // Free memory
        delete[] inMST;
        delete[] key;
        delete[] parent;

        return result;
    }

    /**
     * @brief Kruskal's algorithm for finding Minimum Spanning Tree
     *
     * Implements Kruskal's algorithm to find a Minimum Spanning Tree (MST)
     * of the given graph. Requires at least n-1 edges to form an MST.
     * Uses Union-Find data structure for cycle detection.
     *
     * @param graph The input graph
     * @return Graph A tree representing the MST
     * @throws std::runtime_error If the graph doesn't have at least n-1 edges
     */
    Graph Algorithms::kruskal(const Graph &graph)
    {
        // throw an error if the graph has not has at least n-1 edges
        if (graph.getNumEdges() < graph.getNumVertices() - 1)
        {
            throw std::runtime_error("Graph does not have at least n-1 edges!");
        }

        int numVertices = graph.getNumVertices();

        // Create an empty result graph
        Graph result(numVertices);

        // Collect all edges from the graph
        // In the worst case, there are n(n-1)/2 edges in an undirected graph
        int maxEdges = numVertices * (numVertices - 1) / 2;
        Edge *edges = new Edge[maxEdges];
        int edgeCount = 0;

        // Iterate through all vertices and edges
        for (int i = 0; i < numVertices; i++)
        {
            Node *neighbor = graph.getHead(i);
            while (neighbor != nullptr)
            {
                // Prevent duplicates (since this is an undirected graph)
                // Add the edge only if i < neighbor
                if (i < neighbor->val)
                {
                    edges[edgeCount].src = i;
                    edges[edgeCount].dest = neighbor->val;
                    edges[edgeCount].weight = neighbor->cost;
                    edgeCount++;
                }

                // Move to the next neighbor
                neighbor = neighbor->next;
            }
        }

        // Sort edges by weight in ascending order
        for (int i = 0; i < edgeCount - 1; i++)
        {
            for (int j = 0; j < edgeCount - i - 1; j++)
            {
                if (edges[j].weight > edges[j + 1].weight)
                {
                    // Swap
                    Edge temp = edges[j];
                    edges[j] = edges[j + 1];
                    edges[j + 1] = temp;
                }
            }
        }

        // Create Union-Find data structure
        UnionFind uf(numVertices);

        // Process edges in ascending order of weight
        for (int i = 0; i < edgeCount; i++)
        {
            int src = edges[i].src;
            int dest = edges[i].dest;

            // Check if adding this edge creates a cycle
            if (!uf.connected(src, dest))
            {
                // If not, add it to the MST
                result.addEdge(src, dest, edges[i].weight);

                // Union the sets
                uf.unionSet(src, dest);
            }
        }

        // Free memory
        delete[] edges;

        return result;
    }

    // bool Algorithms::isConnected(const Graph &graph)
    // {
    //     int numVertices = graph.getNumVertices();

    //     // Create an empty result graph
    //     Graph result(numVertices);

    //     // Array to mark visited vertices
    //     bool *visited = new bool[numVertices]();

    //     // Perform DFS from the first vertex
    //     dfsUtil(graph, 0, visited, result, -1);

    //     // Check if all vertices were included in the tree
    //     for (int i = 0; i < numVertices; i++)
    //     {
    //         if (!visited[i])
    //         {
    //             std::cout << "Graph is not connected!" << std::endl;
    //             // Free memory
    //             delete[] visited;
    //             return false;
    //         }
    //     }

    //     // Free memory
    //     delete[] visited;

    //     std::cout << "Graph is connected!" << std::endl;
    //     return true;
    // }
}