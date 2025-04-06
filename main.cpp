// 313610123
#include <iostream>
#include <iomanip>
#include "graph.hpp"

using namespace graph;

/**
 * @brief Prints a formatted header with the given title
 *
 * Creates a visually distinct section header by surrounding the title
 * with equal signs and centering it within a fixed width.
 *
 * @param title The title text to display in the header
 */
void printHeader(const std::string &title)
{
    std::cout << "\n\n"
              << std::string(50, '=') << std::endl;
    std::cout << std::setw(25 + title.length() / 2) << title << std::endl;
    std::cout << std::string(50, '=') << std::endl;
}

/**
 * @brief Creates a predefined example graph for demonstration
 *
 * Initializes a graph with 6 vertices and 9 edges with specific weights.
 * This provides a consistent example for demonstrating the various graph algorithms.
 *
 * @return Graph A predefined graph with 6 vertices and 9 edges
 */
Graph createExampleGraph()
{
    // Create a graph with 6 vertices
    int numVertices = 6;
    int numEdges = 9;

    Edge edges[9] = {
        {0, 1, 4}, // Edge from 0 to 1 with weight 4
        {0, 2, 3}, // Edge from 0 to 2 with weight 3
        {1, 2, 2}, // Edge from 1 to 2 with weight 2
        {1, 3, 5}, // Edge from 1 to 3 with weight 5
        {2, 3, 7}, // Edge from 2 to 3 with weight 7
        {2, 4, 8}, // Edge from 2 to 4 with weight 8
        {3, 4, 1}, // Edge from 3 to 4 with weight 1
        {3, 5, 6}, // Edge from 3 to 5 with weight 6
        {4, 5, 9}  // Edge from 4 to 5 with weight 9
    };

    return Graph(edges, numEdges, numVertices);
}

/**
 * @brief Creates a user-defined graph through console input
 *
 * Prompts the user to enter the number of vertices and edges,
 * and then collects information for each edge (source, destination, weight).
 * Performs validation to ensure the graph is valid (no self-loops or duplicate edges).
 *
 * @return Graph A graph created based on user input
 * @throws std::invalid_argument If the user provides invalid input
 */
Graph createUserGraph()
{
    int numVertices;

    std::cout << "Enter number of vertices: ";
    std::cin >> numVertices;
    if (numVertices <= 0)
    {
        throw std::invalid_argument("Invalid number of vertices!");
    }

    // Create an empty graph with the specified number of vertices
    Graph graph(numVertices);

    std::cout << "Enter number of edges: ";
    int numEdges;
    std::cin >> numEdges;
    if (numEdges <= 0)
    {
        throw std::invalid_argument("Invalid number of edges!");
    }

    std::cout << "\nFor each edge enter: source vertex, destination vertex, weight\n";
    for (int i = 0; i < numEdges; i++)
    {
        int src, dest, weight;
        std::cout << "Edge " << i + 1 << ": ";
        std::cin >> src >> dest >> weight;

        // First perform checks that are not included in addEdge
        // or those for which we want to give the user specific error messages

        // Check validity of vertex numbers
        if (src >= numVertices || dest >= numVertices || src < 0 || dest < 0)
        {
            std::cout << "Invalid vertex numbers! Vertices should be between 0 and "
                      << numVertices - 1 << std::endl;
            i--; // Try again
            continue;
        }

        // Check for self-loops
        if (src == dest)
        {
            std::cout << "Self-loops are not allowed in a simple graph! Please try again." << std::endl;
            i--; // Try again
            continue;
        }

        // Check for duplicate edges
        if (graph.hasEdge(src, dest))
        {
            std::cout << "Duplicate edge detected! A simple graph cannot have multiple edges between the same vertices." << std::endl;
            i--; // Try again
            continue;
        }

        // If we've reached here, the edge is valid and can be added to the graph
        try
        {
            graph.addEdge(src, dest, weight);
        }
        catch (const std::exception &e)
        {
            // In case the addEdge function throws an unexpected exception
            std::cout << "Error adding edge: " << e.what() << ". Please try again." << std::endl;
            i--; // Try again
            continue;
        }
    }

    return graph;
}

/**
 * @brief Main function that demonstrates the graph algorithms
 *
 * Provides a menu for the user to choose between a predefined example graph
 * or creating their own graph. Then runs all the graph algorithms (BFS, DFS,
 * Dijkstra, Prim, and Kruskal) on the selected graph and displays the results.
 *
 * @return int Exit status code (0 for success, 1 for error)
 */
int main()
{
    std::cout << "Graph Algorithms Demonstration" << std::endl;
    std::cout << "----------------------------" << std::endl;

    int choice;
    std::cout << "Choose an option:" << std::endl;
    std::cout << "1. Use example graph" << std::endl;
    std::cout << "2. Create your own graph" << std::endl;
    std::cout << "Your choice: ";
    std::cin >> choice;

    if (choice != 1 && choice != 2)
    {
        std::cerr << "Invalid choice! Please try again." << std::endl;
        return 1;
    }

    Graph graph = (choice == 1) ? createExampleGraph() : createUserGraph();

    // Print the original graph
    printHeader("ORIGINAL GRAPH");
    std::cout << "Graph structure with " << graph.getNumVertices() << " vertices:" << std::endl;
    for (int i = 0; i < graph.getNumVertices(); i++)
    {
        graph.print_graph(graph.getHead(i), i);
    }

    // Choose a starting vertex for algorithms
    int startVertex;
    std::cout << "\nChoose a start vertex for BFS, DFS and Dijkstra (0-" << graph.getNumVertices() - 1 << "): ";
    std::cin >> startVertex;

    // Validate the starting vertex
    if (startVertex < 0 || startVertex >= graph.getNumVertices())
    {
        std::cout << "Invalid vertex! Using vertex 0 instead." << std::endl;
        startVertex = 0;
    }

    try
    {
        // BFS Algorithm
        printHeader("BFS ALGORITHM");
        std::cout << "Running BFS starting from vertex " << startVertex << ":" << std::endl;
        Graph bfsResult = Algorithms::bfs(graph, startVertex);
        std::cout << "BFS result (tree structure):" << std::endl;
        for (int i = 0; i < bfsResult.getNumVertices(); i++)
        {
            bfsResult.print_graph(bfsResult.getHead(i), i);
        }

        // DFS Algorithm
        printHeader("DFS ALGORITHM");
        std::cout << "Running DFS starting from vertex " << startVertex << ":" << std::endl;
        Graph dfsResult = Algorithms::dfs(graph, startVertex);
        std::cout << "DFS result (tree structure):" << std::endl;
        for (int i = 0; i < dfsResult.getNumVertices(); i++)
        {
            dfsResult.print_graph(dfsResult.getHead(i), i);
        }

        // Dijkstra's Algorithm
        printHeader("DIJKSTRA ALGORITHM");
        std::cout << "Running Dijkstra's algorithm from vertex " << startVertex << ":" << std::endl;
        Graph dijkstraResult = Algorithms::dijkstra(graph, startVertex);
        std::cout << "Dijkstra result (shortest paths tree):" << std::endl;
        for (int i = 0; i < dijkstraResult.getNumVertices(); i++)
        {
            dijkstraResult.print_graph(dijkstraResult.getHead(i), i);
        }

        // Prim's Algorithm
        printHeader("PRIM'S ALGORITHM");
        std::cout << "Running Prim's algorithm for Minimum Spanning Tree:" << std::endl;
        Graph primResult = Algorithms::prim(graph);
        std::cout << "Prim result (MST):" << std::endl;
        for (int i = 0; i < primResult.getNumVertices(); i++)
        {
            primResult.print_graph(primResult.getHead(i), i);
        }

        // Kruskal's Algorithm
        printHeader("KRUSKAL'S ALGORITHM");
        std::cout << "Running Kruskal's algorithm for Minimum Spanning Tree:" << std::endl;
        Graph kruskalResult = Algorithms::kruskal(graph);
        std::cout << "Kruskal result (MST):" << std::endl;
        for (int i = 0; i < kruskalResult.getNumVertices(); i++)
        {
            kruskalResult.print_graph(kruskalResult.getHead(i), i);
        }

        // // Edge operations demonstration
        // printHeader("EDGE OPERATIONS");
        // std::cout << "Demonstrating edge operations:" << std::endl;

        // // Copy the original graph
        // Graph testGraph = createExampleGraph();

        // // Print the original graph
        // std::cout << "Original graph:" << std::endl;
        // for (int i = 0; i < testGraph.getNumVertices(); i++)
        // {
        //     testGraph.print_graph(testGraph.getHead(i), i);
        // }

        // // Add a new edge
        // std::cout << "\nAdding edge from 0 to 5 with weight 10:" << std::endl;
        // testGraph.addEdge(0, 5, 10);
        // std::cout << "After adding edge:" << std::endl;
        // for (int i = 0; i < testGraph.getNumVertices(); i++)
        // {
        //     testGraph.print_graph(testGraph.getHead(i), i);
        // }

        // // Remove an edge
        // std::cout << "\nRemoving edge between 0 and 1:" << std::endl;
        // testGraph.removeEdge(0, 1);
        // std::cout << "After removing edge:" << std::endl;
        // for (int i = 0; i < testGraph.getNumVertices(); i++)
        // {
        //     testGraph.print_graph(testGraph.getHead(i), i);
        // }

        // // Check if edges exist
        // std::cout << "\nChecking if edges exist:" << std::endl;
        // std::cout << "Edge 0-2 exists: " << (testGraph.hasEdge(0, 2) ? "Yes" : "No") << std::endl;
        // std::cout << "Edge 0-1 exists: " << (testGraph.hasEdge(0, 1) ? "Yes" : "No") << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    std::cout << "\nDemonstration completed successfully!" << std::endl;
    return 0;
}