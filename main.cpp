#include <iostream>
#include "graph.hpp"

using namespace graph;

int main()
{
    int numVertices;

    std::cout << "Enter number of vertices: ";
    std::cin >> numVertices;

    std::cout << "Enter number of edges: ";
    int numEdges;
    std::cin >> numEdges;

    Edge *edges = new Edge[numEdges];
    int validEdges = 0; // Track number of valid edges

    std::cout << "\nFor each edge enter: source vertex, destination vertex, weight\n";
    for (int i = 0; i < numEdges; i++)
    {
        std::cout << "Edge " << i + 1 << ": ";
        std::cin >> edges[i].src >> edges[i].dest >> edges[i].weight;

        // Validate input - check vertex numbers
        if (edges[i].src >= numVertices || edges[i].dest >= numVertices ||
            edges[i].src < 0 || edges[i].dest < 0)
        {
            std::cout << "Invalid vertex numbers! Vertices should be between 0 and "
                      << numVertices - 1 << std::endl;
            i--; // Try again with same array position
            continue;
        }

        // Check for self-loops
        if (edges[i].src == edges[i].dest)
        {
            std::cout << "Self-loops are not allowed in a simple graph! Please try again." << std::endl;
            i--; // Try again with same array position
            continue;
        }

        // Check for duplicate edges
        bool isDuplicate = false;
        for (int j = 0; j < validEdges; j++)
        {
            if ((edges[j].src == edges[i].src && edges[j].dest == edges[i].dest) ||
                (edges[j].src == edges[i].dest && edges[j].dest == edges[i].src))
            {
                std::cout << "Duplicate edge detected! A simple graph cannot have multiple edges between the same vertices." << std::endl;
                isDuplicate = true;
                break;
            }
        }

        if (isDuplicate)
        {
            i--; // Try again with same array position
            continue;
        }

        // If we got here, the edge is valid
        validEdges++;
    }

    // Create the graph
    Graph graph(edges, validEdges, numVertices);

    // Print the graph
    std::cout << "\nGraph structure:\n";
    for (int i = 0; i < numVertices; i++)
    {
        graph.print_graph(graph.head[i], i);
    }

    delete[] edges;
    return 0;
}