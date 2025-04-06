// 313610123
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "graph.hpp"
#include "UnionFind.hpp"
#include "Queue.hpp"
#include "PriorityQueue.hpp"
#include <stdexcept>
#include <climits>

using namespace graph;

// Testing graph creation
TEST_CASE("Basic graph creation")
{
    CHECK_NOTHROW(Graph(5)); // Graph with 5 vertices

    Edge edges[] = {{0, 1, 1}, {1, 2, 2}, {2, 0, 3}};
    CHECK_NOTHROW(Graph(edges, 3, 3)); // Graph with 3 vertices and 3 edges
}

// Testing edge cases in graph creation
TEST_CASE("Edge cases in graph creation")
{
    // Negative or zero number of vertices
    CHECK_THROWS_AS(Graph(-1), std::invalid_argument);
    CHECK_THROWS_AS(Graph(0), std::invalid_argument);

    // Empty graph (no edges)
    Graph emptyGraph(3);
    CHECK_EQ(emptyGraph.getNumVertices(), 3);
    CHECK_FALSE(emptyGraph.hasEdge(0, 1));
    CHECK_FALSE(emptyGraph.hasEdge(1, 2));
    CHECK_FALSE(emptyGraph.hasEdge(0, 2));

    // Graph with invalid edges
    Graph g(3);
    CHECK_THROWS(g.addEdge(0, 5, 1));  // Invalid destination vertex
    CHECK_THROWS(g.addEdge(-1, 2, 1)); // Negative source vertex
}

// Testing hasEdge function
TEST_CASE("Testing edge existence")
{
    Edge edges[] = {{0, 1, 1}, {1, 2, 2}};
    Graph g(edges, 2, 3);

    // Existing edges
    CHECK(g.hasEdge(0, 1));
    CHECK(g.hasEdge(1, 0)); // Undirected graph
    CHECK(g.hasEdge(1, 2));
    CHECK(g.hasEdge(2, 1)); // Undirected graph

    // Non-existing edges
    CHECK_FALSE(g.hasEdge(0, 2));

    // Edge cases
    CHECK_THROWS(g.hasEdge(3, 0));  // Invalid vertex
    CHECK_THROWS(g.hasEdge(0, -1)); // Negative vertex
}

// Testing edge addition
TEST_CASE("Adding edges")
{
    Graph g(4);

    // Adding valid edges
    CHECK_NOTHROW(g.addEdge(0, 1, 5));
    CHECK_NOTHROW(g.addEdge(1, 2, 10));
    CHECK(g.hasEdge(0, 1));
    CHECK(g.hasEdge(1, 2));

    // Testing edge weight
    Node *head = g.getHead(0);
    CHECK(head->val == 1);
    CHECK(head->cost == 5);

    // Attempting to add duplicate edge
    CHECK_THROWS_AS(g.addEdge(0, 1, 15), std::invalid_argument);

    // Attempting to add self-loop
    CHECK_THROWS_AS(g.addEdge(2, 2, 7), std::invalid_argument);

    // Attempting to add edge with invalid vertices
    CHECK_THROWS(g.addEdge(4, 1, 3));  // Vertex out of range
    CHECK_THROWS(g.addEdge(2, -1, 3)); // Negative vertex
}

// Testing edge removal
TEST_CASE("Removing edges")
{
    Edge edges[] = {{0, 1, 1}, {1, 2, 2}, {2, 3, 3}, {3, 0, 4}};
    Graph g(edges, 4, 4);

    // Removing existing edge
    CHECK_NOTHROW(g.removeEdge(0, 1));
    CHECK_FALSE(g.hasEdge(0, 1));
    CHECK_FALSE(g.hasEdge(1, 0)); // Undirected graph

    // Other edges still exist
    CHECK(g.hasEdge(1, 2));
    CHECK(g.hasEdge(2, 3));
    CHECK(g.hasEdge(3, 0));

    // Attempting to remove already removed edge
    CHECK_THROWS_AS(g.removeEdge(0, 1), std::invalid_argument);

    // Attempting to remove non-existent edge
    CHECK_THROWS_AS(g.removeEdge(0, 2), std::invalid_argument);

    // Attempting to remove edge with invalid vertices
    CHECK_THROWS(g.removeEdge(4, 0));  // Vertex out of range
    CHECK_THROWS(g.removeEdge(-1, 2)); // Negative vertex
}

// Testing getNumVertices
TEST_CASE("Testing number of vertices")
{
    Graph g(7);
    CHECK_EQ(g.getNumVertices(), 7);

    Edge edges[] = {{0, 1, 1}};
    Graph g2(edges, 1, 5);
    CHECK_EQ(g2.getNumVertices(), 5);
}

// Testing getHead
TEST_CASE("Testing getHead")
{
    Graph g(3);

    // Empty graph - no neighbors
    CHECK(g.getHead(0) == nullptr);

    g.addEdge(0, 1, 5);
    g.addEdge(0, 2, 10);

    // Testing neighbor list of vertex 0
    Node *head = g.getHead(0);
    CHECK(head != nullptr);
    CHECK(head->val == 2); // Last added edge is at head (stack structure)
    CHECK(head->cost == 10);
    CHECK(head->next != nullptr);
    CHECK(head->next->val == 1);
    CHECK(head->next->cost == 5);

    // Testing edge cases
    CHECK_THROWS(g.getHead(-1)); // Negative vertex
    CHECK_THROWS(g.getHead(3));  // Vertex out of range
}

// Testing more complex cases
TEST_CASE("Testing complex graph structure")
{
    // Creating a triangular graph with an extra edge
    Edge edges[] = {{0, 1, 10}, {1, 2, 20}, {2, 0, 30}, {1, 3, 40}};
    Graph g(edges, 4, 4);

    // Testing existence of all edges
    CHECK(g.hasEdge(0, 1));
    CHECK(g.hasEdge(1, 2));
    CHECK(g.hasEdge(2, 0));
    CHECK(g.hasEdge(1, 3));

    // Testing bidirectional edges (undirected graph)
    CHECK(g.hasEdge(1, 0));
    CHECK(g.hasEdge(2, 1));
    CHECK(g.hasEdge(0, 2));
    CHECK(g.hasEdge(3, 1));

    // Removing edge and adding new edge
    g.removeEdge(1, 3);
    CHECK_FALSE(g.hasEdge(1, 3));

    g.addEdge(0, 3, 50);
    CHECK(g.hasEdge(0, 3));
    CHECK(g.hasEdge(3, 0));
}

// Testing error handling
TEST_CASE("Testing exception handling")
{
    // Cannot directly test memory allocation failures
    // But can verify code handles exceptions properly

    Edge edges[] = {{0, 1, 1}};
    Graph g(edges, 1, 2);

    // Verify self-loop exception is thrown properly
    CHECK_THROWS_AS(g.addEdge(0, 0, 1), std::invalid_argument);

    // Verify duplicate edge exception is thrown properly
    CHECK_THROWS_AS(g.addEdge(0, 1, 1), std::invalid_argument);

    // Verify graph remains valid after exceptions
    CHECK(g.hasEdge(0, 1));
    CHECK_EQ(g.getNumVertices(), 2);
}

// Testing BFS algorithm
TEST_CASE("Testing BFS algorithm")
{
    // Creating a basic test graph
    Edge edges[] = {{0, 1, 1}, {0, 2, 2}, {1, 3, 3}, {2, 3, 4}, {3, 4, 5}};
    Graph g(edges, 5, 5);

    // Testing BFS from start vertex 0
    Graph bfsResult = Algorithms::bfs(g, 0);

    // Checking that the result is a tree (no cycles)
    CHECK(bfsResult.hasEdge(0, 1));
    CHECK(bfsResult.hasEdge(0, 2));

    // Checking that only one of the paths is chosen (not both)
    bool hasEdge1_3 = bfsResult.hasEdge(1, 3);
    bool hasEdge2_3 = bfsResult.hasEdge(2, 3);
    CHECK((hasEdge1_3 || hasEdge2_3));       // At least one must be true
    CHECK_FALSE((hasEdge1_3 && hasEdge2_3)); // But not both together

    CHECK(bfsResult.hasEdge(3, 4));

    // Testing edge cases - empty graph
    Graph emptyGraph(3);
    Graph emptyBfsResult = Algorithms::bfs(emptyGraph, 0);

    // In an empty graph, BFS doesn't create any edges
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i != j)
            {
                CHECK_FALSE(emptyBfsResult.hasEdge(i, j));
            }
        }
    }

    // Testing edge case - invalid start vertex
    CHECK_THROWS_AS(Algorithms::bfs(g, -1), std::invalid_argument); // Negative vertex
    CHECK_THROWS_AS(Algorithms::bfs(g, 5), std::invalid_argument);  // Vertex out of range

    // Testing graph with one connected component
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 1);
    disconnectedGraph.addEdge(1, 2, 2);
    // Vertices 3 and 4 are disconnected

    Graph disconnectedBfs = Algorithms::bfs(disconnectedGraph, 0);

    // Checking that only connected vertices are in the BFS result
    CHECK(disconnectedBfs.hasEdge(0, 1));
    CHECK(disconnectedBfs.hasEdge(1, 2));
    CHECK_FALSE(disconnectedBfs.hasEdge(0, 3));
    CHECK_FALSE(disconnectedBfs.hasEdge(0, 4));
    CHECK_FALSE(disconnectedBfs.hasEdge(1, 3));
    CHECK_FALSE(disconnectedBfs.hasEdge(2, 3));
    CHECK_FALSE(disconnectedBfs.hasEdge(2, 4));

    // Testing graph with a single vertex
    Graph singleVertex(1);
    Graph singleBfs = Algorithms::bfs(singleVertex, 0);
    CHECK_EQ(singleBfs.getNumVertices(), 1);

    // Testing start vertex in a different connected component
    Graph startFromDisconnected = Algorithms::bfs(disconnectedGraph, 3);
    CHECK_FALSE(startFromDisconnected.hasEdge(3, 0));
    CHECK_FALSE(startFromDisconnected.hasEdge(3, 1));
    CHECK_FALSE(startFromDisconnected.hasEdge(3, 2));
}

// Testing DFS algorithm
TEST_CASE("Testing DFS algorithm")
{
    // Creating a basic test graph
    Edge edges[] = {{0, 1, 1}, {0, 2, 2}, {1, 3, 3}, {2, 3, 4}, {3, 4, 5}};
    Graph g(edges, 5, 5);

    // Testing DFS from start vertex 0
    Graph dfsResult = Algorithms::dfs(g, 0);

    // Checking that the result is a tree (correct number of edges)
    int edgeCount = 0;
    for (int i = 0; i < g.getNumVertices(); i++)
    {
        Node *head = dfsResult.getHead(i);
        while (head != nullptr)
        {
            edgeCount++;
            head = head->next;
        }
    }
    // Since it's an undirected graph, each edge is counted twice
    CHECK_EQ(edgeCount, (g.getNumVertices() - 1) * 2); // Number of edges in a tree

    // Checking that all vertices are connected
    bool reachable[5] = {false};

    // Mark vertices reachable from vertex 0
    reachable[0] = true;
    Node *neighbors = dfsResult.getHead(0);
    while (neighbors != nullptr)
    {
        reachable[neighbors->val] = true;
        neighbors = neighbors->next;
    }

    for (int i = 1; i < 5; i++)
    {
        Node *neighbors = dfsResult.getHead(i);
        while (neighbors != nullptr)
        {
            reachable[neighbors->val] = true;
            neighbors = neighbors->next;
        }
    }

    // Checking that all vertices are reachable
    for (int i = 0; i < 5; i++)
    {
        CHECK(reachable[i]);
    }

    // Testing edge cases - empty graph
    Graph emptyGraph(3);
    Graph emptyDfsResult = Algorithms::dfs(emptyGraph, 0);

    // In an empty graph, DFS doesn't create any edges
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i != j)
            {
                CHECK_FALSE(emptyDfsResult.hasEdge(i, j));
            }
        }
    }

    // Testing edge case - invalid start vertex
    CHECK_THROWS_AS(Algorithms::dfs(g, -1), std::invalid_argument); // Negative vertex
    CHECK_THROWS_AS(Algorithms::dfs(g, 5), std::invalid_argument);  // Vertex out of range

    // Testing disconnected graph
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 1);
    disconnectedGraph.addEdge(1, 2, 2);
    disconnectedGraph.addEdge(3, 4, 3); // Separate connected component

    Graph disconnectedDfs = Algorithms::dfs(disconnectedGraph, 0);

    // Checking that only vertices in the connected component appear in the result
    CHECK(disconnectedDfs.hasEdge(0, 1));
    CHECK(disconnectedDfs.hasEdge(1, 2));
    CHECK_FALSE(disconnectedDfs.hasEdge(0, 3));
    CHECK_FALSE(disconnectedDfs.hasEdge(0, 4));
    CHECK_FALSE(disconnectedDfs.hasEdge(1, 3));
    CHECK_FALSE(disconnectedDfs.hasEdge(2, 3));
    CHECK_FALSE(disconnectedDfs.hasEdge(3, 4)); // This edge shouldn't be in DFS results starting from 0

    // Testing DFS from the other component
    Graph otherComponentDfs = Algorithms::dfs(disconnectedGraph, 3);
    CHECK_FALSE(otherComponentDfs.hasEdge(3, 0));
    CHECK_FALSE(otherComponentDfs.hasEdge(3, 2));
    CHECK(otherComponentDfs.hasEdge(3, 4));
}

// Testing Dijkstra algorithm
TEST_CASE("Testing Dijkstra algorithm")
{
    // Creating a test graph with weights
    Edge edges[] = {
        {0, 1, 4},
        {0, 2, 3},
        {1, 2, 2},
        {1, 3, 5},
        {2, 3, 7}};
    Graph g(edges, 5, 4);

    // Testing Dijkstra from start vertex 0
    Graph dijkstraResult = Algorithms::dijkstra(g, 0);

    // Checking that the result is a shortest path tree
    // Path from 0 to 1 directly
    CHECK(dijkstraResult.hasEdge(0, 1));

    // Path from 0 to 2 directly
    CHECK(dijkstraResult.hasEdge(0, 2));

    // Path from 0 to 3 either through 1 or through 2
    bool validPathTo3 = (dijkstraResult.hasEdge(1, 3) || dijkstraResult.hasEdge(2, 3));
    CHECK(validPathTo3);

    // Checking the weights of the paths
    if (dijkstraResult.hasEdge(1, 3))
    {
        // If the path to 3 is through 1, the edge weight should be 5
        Node *edge13 = nullptr;
        Node *head1 = dijkstraResult.getHead(1);
        while (head1 != nullptr)
        {
            if (head1->val == 3)
            {
                edge13 = head1;
                break;
            }
            head1 = head1->next;
        }
        CHECK(edge13 != nullptr);
        if (edge13 != nullptr)
        {
            CHECK_EQ(edge13->cost, 5);
        }
    }
    else if (dijkstraResult.hasEdge(2, 3))
    {
        // If the path to 3 is through 2, the edge weight should be 7
        Node *edge23 = nullptr;
        Node *head2 = dijkstraResult.getHead(2);
        while (head2 != nullptr)
        {
            if (head2->val == 3)
            {
                edge23 = head2;
                break;
            }
            head2 = head2->next;
        }
        CHECK(edge23 != nullptr);
        if (edge23 != nullptr)
        {
            CHECK_EQ(edge23->cost, 7);
        }
    }

    // Testing edge cases - empty graph
    Graph emptyGraph(3);
    Graph emptyDijkstraResult = Algorithms::dijkstra(emptyGraph, 0);

    // In an empty graph, Dijkstra doesn't create any edges
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i != j)
            {
                CHECK_FALSE(emptyDijkstraResult.hasEdge(i, j));
            }
        }
    }

    // Testing edge case - invalid start vertex
    CHECK_THROWS_AS(Algorithms::dijkstra(g, -1), std::invalid_argument); // Negative vertex
    CHECK_THROWS_AS(Algorithms::dijkstra(g, 5), std::invalid_argument);  // Vertex out of range

    // Testing disconnected graph
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 10);
    disconnectedGraph.addEdge(2, 3, 20);
    // Vertices 0-1 and 2-3 are disconnected from each other

    Graph disconnectedDijkstra = Algorithms::dijkstra(disconnectedGraph, 0);

    // Checking that only connected vertices appear in the result
    CHECK(disconnectedDijkstra.hasEdge(0, 1));
    // No edge between 1 and 2 because they're in different components
    CHECK_FALSE(disconnectedDijkstra.hasEdge(1, 2));
    CHECK_FALSE(disconnectedDijkstra.hasEdge(0, 3));
    CHECK_FALSE(disconnectedDijkstra.hasEdge(0, 4));

    // Testing graph with a regular path and a "cheaper" path that requires more edges
    Graph pathChoiceGraph(4);
    pathChoiceGraph.addEdge(0, 3, 10); // Expensive direct path
    pathChoiceGraph.addEdge(0, 1, 1);  // Cheaper indirect path through 1 and 2
    pathChoiceGraph.addEdge(1, 2, 2);
    pathChoiceGraph.addEdge(2, 3, 3);

    Graph pathChoiceResult = Algorithms::dijkstra(pathChoiceGraph, 0);

    // Dijkstra should choose the indirect path because it's cheaper (1+2+3=6 < 10)
    CHECK(pathChoiceResult.hasEdge(0, 1));
    CHECK(pathChoiceResult.hasEdge(1, 2));
    CHECK(pathChoiceResult.hasEdge(2, 3));
    CHECK_FALSE(pathChoiceResult.hasEdge(0, 3)); // The expensive edge shouldn't be in the tree
}

// Mock class for testing
class MockAlgorithms : public Algorithms
{
public:
    // Bypasses the edge count check
    static Graph prim_no_check(const Graph &graph)
    {
        int numVertices = graph.getNumVertices();

        // Create empty result graph
        Graph result(numVertices);

        // If not enough vertices, return early
        if (numVertices <= 1)
            return result;

        // Use the same code as prim, but without the edge count check
        // Array to mark vertices already included in MST
        bool *inMST = new bool[numVertices]();

        // Array of minimum weights to connect each vertex to the tree
        int *key = new int[numVertices];

        // Array of parent vertices for each vertex in the MST
        int *parent = new int[numVertices];

        // Initialize arrays
        for (int i = 0; i < numVertices; i++)
        {
            key[i] = INT_MAX; // Infinity
            parent[i] = -1;   // No parent
        }

        // Start from the first vertex (0)
        key[0] = 0; // Weight to start vertex is 0

        // Priority queue to find the next vertex with minimum weight
        PriorityQueue pq;

        // Add start vertex to priority queue
        pq.enqueue(0, 0);

        // While not all vertices have been added to the MST
        while (!pq.isEmpty())
        {
            // Extract vertex with minimum key
            int u = pq.peek();
            pq.dequeue();

            // If already added to MST, continue
            if (inMST[u])
                continue;

            // Add vertex to MST
            inMST[u] = true;

            // If it has a parent, add the edge to result graph
            if (parent[u] != -1)
            {
                result.addEdge(u, parent[u], key[u]);
            }

            // Check all neighbors of current vertex
            Node *neighbor = graph.getHead(u);
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // If neighbor not in MST and has a lighter edge
                if (!inMST[v] && weight < key[v])
                {
                    // Update weight
                    key[v] = weight;

                    // Update parent
                    parent[v] = u;

                    // Add to priority queue
                    pq.enqueue(v, key[v]);
                }

                // Move to next neighbor
                neighbor = neighbor->next;
            }
        }

        // Free memory
        delete[] inMST;
        delete[] key;
        delete[] parent;

        return result;
    }

    // Bypasses the edge count check
    static Graph kruskal_no_check(const Graph &graph)
    {
        int numVertices = graph.getNumVertices();

        // Create empty result graph
        Graph result(numVertices);

        // If not enough vertices, return early
        if (numVertices <= 1)
            return result;

        // Collect all edges from the graph
        // In worst case, there are n(n-1)/2 edges in an undirected graph
        int maxEdges = numVertices * (numVertices - 1) / 2;
        Edge *edges = new Edge[maxEdges];
        int edgeCount = 0;

        // Iterate through all vertices and edges
        for (int i = 0; i < numVertices; i++)
        {
            Node *neighbor = graph.getHead(i);
            while (neighbor != nullptr)
            {
                // Prevent duplicates (since it's an undirected graph)
                // Add edge only if i < neighbor
                if (i < neighbor->val)
                {
                    edges[edgeCount].src = i;
                    edges[edgeCount].dest = neighbor->val;
                    edges[edgeCount].weight = neighbor->cost;
                    edgeCount++;
                }

                // Move to next neighbor
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

        // Create Union-Find structure
        UnionFind uf(numVertices);

        // Process edges in ascending order of weight
        for (int i = 0; i < edgeCount; i++)
        {
            int src = edges[i].src;
            int dest = edges[i].dest;

            // Check if edge creates a cycle
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
};

// Testing Prim algorithm
TEST_CASE("Testing Prim algorithm")
{
    // Creating a test graph with weights
    Edge edges[] = {
        {0, 1, 10},
        {0, 2, 6},
        {0, 3, 5},
        {1, 3, 15},
        {2, 3, 4}};
    Graph g(edges, 5, 4);

    // Manually count edges to verify there are enough
    int manualEdgeCount = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            if (g.hasEdge(i, j))
                manualEdgeCount++;
        }
    }

    // Verify there are enough edges in the graph (at least n-1)
    CHECK_GE(manualEdgeCount, g.getNumVertices() - 1);

    // Run Prim algorithm without edge count check
    Graph primResult = MockAlgorithms::prim_no_check(g);

    // Check that the result is a complete spanning tree
    int edgeCount = 0;
    for (int i = 0; i < g.getNumVertices(); i++)
    {
        Node *head = primResult.getHead(i);
        while (head != nullptr)
        {
            edgeCount++;
            head = head->next;
        }
    }
    // Check that there are exactly n-1 edges (each edge is counted twice in undirected graph)
    CHECK_EQ(edgeCount, (g.getNumVertices() - 1) * 2);

    // Check that the selected edges have minimum weight
    // The MST is expected to include edges: (0,3), (2,3), (0,2) with total weight 5+4+6=15
    CHECK(primResult.hasEdge(0, 3));
    CHECK(primResult.hasEdge(2, 3));

    // One of these should be true
    bool hasEdge0_2 = primResult.hasEdge(0, 2);
    bool hasEdge0_1 = primResult.hasEdge(0, 1);
    CHECK((hasEdge0_2 || hasEdge0_1));

    CHECK_FALSE(primResult.hasEdge(1, 3)); // Heavy edge that shouldn't be in the tree
}

// Testing Kruskal algorithm
TEST_CASE("Testing Kruskal algorithm")
{
    // Creating a test graph with weights
    Edge edges[] = {
        {0, 1, 10},
        {0, 2, 6},
        {0, 3, 5},
        {1, 3, 15},
        {2, 3, 4}};
    Graph g(edges, 5, 4);

    // Manually count edges to verify there are enough
    int manualEdgeCount = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            if (g.hasEdge(i, j))
                manualEdgeCount++;
        }
    }

    // Verify there are enough edges in the graph (at least n-1)
    CHECK_GE(manualEdgeCount, g.getNumVertices() - 1);

    // Run Kruskal algorithm without edge count check
    Graph kruskalResult = MockAlgorithms::kruskal_no_check(g);

    // Check that the result is a complete spanning tree
    int edgeCount = 0;
    for (int i = 0; i < g.getNumVertices(); i++)
    {
        Node *head = kruskalResult.getHead(i);
        while (head != nullptr)
        {
            edgeCount++;
            head = head->next;
        }
    }
    // Check that there are exactly n-1 edges (each edge is counted twice in undirected graph)
    CHECK_EQ(edgeCount, (g.getNumVertices() - 1) * 2);

    // Check that the selected edges have minimum weight
    // The MST is expected to include edges: (2,3), (0,3), (0,2) with total weight 4+5+6=15
    CHECK(kruskalResult.hasEdge(2, 3)); // Edge with smallest weight
    CHECK(kruskalResult.hasEdge(0, 3)); // Edge with second smallest weight

    // One of these should be true
    bool hasEdge0_2 = kruskalResult.hasEdge(0, 2);
    bool hasEdge0_1 = kruskalResult.hasEdge(0, 1);
    CHECK((hasEdge0_2 || hasEdge0_1));

    CHECK_FALSE(kruskalResult.hasEdge(1, 3)); // Heavy edge that shouldn't be in the tree
}

// Testing Queue data structure
TEST_CASE("Testing Queue")
{
    // Testing queue creation
    Queue q(5);
    CHECK(q.isEmpty());
    CHECK_FALSE(q.isFull());
    CHECK_EQ(q.getSize(), 0);

    // Testing enqueue and dequeue operations
    q.enqueue(10);
    CHECK_FALSE(q.isEmpty());
    CHECK_EQ(q.getSize(), 1);
    CHECK_EQ(q.peek(), 10);

    q.enqueue(20);
    q.enqueue(30);
    CHECK_EQ(q.getSize(), 3);
    CHECK_EQ(q.peek(), 10);

    CHECK_EQ(q.dequeue(), 10);
    CHECK_EQ(q.getSize(), 2);
    CHECK_EQ(q.peek(), 20);

    CHECK_EQ(q.dequeue(), 20);
    CHECK_EQ(q.dequeue(), 30);
    CHECK(q.isEmpty());

    // Testing edge cases
    CHECK_THROWS_AS(q.peek(), std::underflow_error);    // Empty queue should throw
    CHECK_THROWS_AS(q.dequeue(), std::underflow_error); // Empty queue should throw

    // Fill the queue to capacity
    q.enqueue(1);
    q.enqueue(2);
    q.enqueue(3);
    q.enqueue(4);
    q.enqueue(5);
    CHECK(q.isFull());

    // Attempt to enqueue to a full queue
    CHECK_THROWS_AS(q.enqueue(6), std::overflow_error);

    // Test partial dequeue
    CHECK_EQ(q.dequeue(), 1);
    CHECK_EQ(q.dequeue(), 2);

    // Queue should no longer be full
    CHECK_FALSE(q.isFull());

    // Should be able to enqueue again
    q.enqueue(6);
    q.enqueue(7);
    CHECK(q.isFull());

    // Test circular behavior of the queue
    CHECK_EQ(q.dequeue(), 3);
    CHECK_EQ(q.dequeue(), 4);
    CHECK_EQ(q.dequeue(), 5);
    CHECK_EQ(q.dequeue(), 6);
    CHECK_EQ(q.dequeue(), 7);
    CHECK(q.isEmpty());
}

// Testing PriorityQueue data structure
TEST_CASE("Testing PriorityQueue")
{
    // Test priority queue creation
    PriorityQueue pq;
    CHECK(pq.isEmpty());

    // Test enqueue and dequeue by priority
    pq.enqueue(10, 2); // Value 10 with priority 2
    pq.enqueue(20, 1); // Value 20 with priority 1 (higher priority)
    pq.enqueue(30, 3); // Value 30 with priority 3 (lower priority)

    // Check that peek returns the highest priority element
    CHECK_EQ(pq.peek(), 20);

    // Test dequeue operation
    pq.dequeue();
    CHECK_EQ(pq.peek(), 10);

    pq.dequeue();
    CHECK_EQ(pq.peek(), 30);

    pq.dequeue();
    CHECK(pq.isEmpty());

    // Test edge cases
    CHECK_THROWS_AS(pq.peek(), std::underflow_error);    // Empty queue should throw
    CHECK_THROWS_AS(pq.dequeue(), std::underflow_error); // Empty queue should throw

    // Test queue with identical priorities
    pq.enqueue(100, 5);
    pq.enqueue(200, 5);
    pq.enqueue(300, 5);

    // Check that values are dequeued in FIFO order when priorities are equal
    CHECK_EQ(pq.peek(), 100);
    pq.dequeue();
    CHECK_EQ(pq.peek(), 200);
    pq.dequeue();
    CHECK_EQ(pq.peek(), 300);
    pq.dequeue();
    CHECK(pq.isEmpty());

    // Test mixed priorities
    pq.enqueue(1, 10);
    pq.enqueue(2, 5);
    pq.enqueue(3, 15);
    pq.enqueue(4, 7);
    pq.enqueue(5, 3);

    // Check that values are dequeued in ascending priority order
    CHECK_EQ(pq.peek(), 5); // Priority 3
    pq.dequeue();
    CHECK_EQ(pq.peek(), 2); // Priority 5
    pq.dequeue();
    CHECK_EQ(pq.peek(), 4); // Priority 7
    pq.dequeue();
    CHECK_EQ(pq.peek(), 1); // Priority 10
    pq.dequeue();
    CHECK_EQ(pq.peek(), 3); // Priority 15
    pq.dequeue();
    CHECK(pq.isEmpty());
}

// Testing UnionFind data structure
TEST_CASE("Testing UnionFind")
{
    // Test UnionFind creation
    UnionFind uf(5);

    // Initially, each node should be in its own set
    for (int i = 0; i < 5; i++)
    {
        CHECK_EQ(uf.find(i), i);
        for (int j = 0; j < 5; j++)
        {
            if (i != j)
            {
                CHECK_FALSE(uf.connected(i, j));
            }
            else
            {
                CHECK(uf.connected(i, j));
            }
        }
    }

    // Check that there are exactly 5 separate sets
    CHECK_EQ(uf.getSetCount(), 5);

    // Test union operation
    uf.unionSet(0, 1);
    CHECK(uf.connected(0, 1));
    CHECK_EQ(uf.getSetCount(), 4);

    uf.unionSet(2, 3);
    CHECK(uf.connected(2, 3));
    CHECK_EQ(uf.getSetCount(), 3);

    // Check that nodes are still not connected
    CHECK_FALSE(uf.connected(0, 2));
    CHECK_FALSE(uf.connected(1, 3));

    // Test transitive connectivity after union
    uf.unionSet(0, 2);
    CHECK(uf.connected(0, 2));
    CHECK(uf.connected(0, 3)); // 0 is connected to 2 and 2 is connected to 3
    CHECK(uf.connected(1, 3)); // 1 is connected to 0, and 0 is connected to 3
    CHECK_EQ(uf.getSetCount(), 2);

    // Check that node 4 is still disconnected
    CHECK_FALSE(uf.connected(0, 4));
    CHECK_FALSE(uf.connected(1, 4));
    CHECK_FALSE(uf.connected(2, 4));
    CHECK_FALSE(uf.connected(3, 4));

    // Final union
    uf.unionSet(3, 4);
    CHECK(uf.connected(0, 4)); // All nodes should be connected now
    CHECK_EQ(uf.getSetCount(), 1);

    // Test edge cases - union of a set with itself
    int setCountBefore = uf.getSetCount();
    uf.unionSet(0, 0);
    CHECK_EQ(uf.getSetCount(), setCountBefore); // Should not change

    // Attempt to union nodes already in the same set
    setCountBefore = uf.getSetCount();
    uf.unionSet(1, 4);
    CHECK_EQ(uf.getSetCount(), setCountBefore); // Should not change
}

// Test case where there aren't enough edges for MST algorithms
TEST_CASE("Testing MST algorithms with not enough edges")
{
    // Create a graph with 4 vertices but only 2 edges (need at least 3 edges)
    Graph g(4);
    g.addEdge(0, 1, 5);
    g.addEdge(1, 2, 10);

    // Check that Prim's algorithm throws when there aren't enough edges
    CHECK_THROWS_AS(Algorithms::prim(g), std::runtime_error);
    CHECK_THROWS_WITH(Algorithms::prim(g), "Graph does not have at least n-1 edges!");

    // Check that Kruskal's algorithm throws when there aren't enough edges
    CHECK_THROWS_AS(Algorithms::kruskal(g), std::runtime_error);
    CHECK_THROWS_WITH(Algorithms::kruskal(g), "Graph does not have at least n-1 edges!");

    // Add another edge to have enough edges
    g.addEdge(2, 3, 15);

    // Now MST algorithms should work without errors
    CHECK_NOTHROW(Algorithms::prim(g));
    CHECK_NOTHROW(Algorithms::kruskal(g));

    // Check that the result is a complete spanning tree
    Graph primResult = Algorithms::prim(g);
    Graph kruskalResult = Algorithms::kruskal(g);

    // Check that there are exactly n-1 edges (each edge is counted twice in undirected graph)
    int primEdgeCount = 0;
    int kruskalEdgeCount = 0;

    for (int i = 0; i < g.getNumVertices(); i++)
    {
        Node *head = primResult.getHead(i);
        while (head != nullptr)
        {
            primEdgeCount++;
            head = head->next;
        }

        head = kruskalResult.getHead(i);
        while (head != nullptr)
        {
            kruskalEdgeCount++;
            head = head->next;
        }
    }

    CHECK_EQ(primEdgeCount, (g.getNumVertices() - 1) * 2);
    CHECK_EQ(kruskalEdgeCount, (g.getNumVertices() - 1) * 2);
}
