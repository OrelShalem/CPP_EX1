#include "graph.hpp"
#include <iostream>

using namespace graph;

Graph::Graph(Edge edges[], int n, int numVertices)
{
    this->head = new Node *[numVertices]();
    this->numVertices = numVertices;

    for (int i = 0; i < numVertices; i++)
    {
        head[i] = nullptr;
    }

    for (int i = 0; i < n; i++)
    {
        int src = edges[i].src;
        int dest = edges[i].dest;
        int weight = edges[i].weight;

        // Check for self-loop (vertex connected to itself)
        if (src == dest)
        {
            std::cout << "Warning: Self-loop detected at vertex " << src << ". Ignoring." << std::endl;
            continue;
        }

        // Check for duplicate edge
        if (hasEdge(src, dest))
        {
            std::cout << "Warning: Duplicate edge detected between " << src << " and "
                      << dest << ". Ignoring." << std::endl;
            continue;
        }

        // add an edge from src to dest
        Node *newNode = getAdjacencyList(dest, weight, head[src]);
        head[src] = newNode;

        // undirected graph
        newNode = getAdjacencyList(src, weight, head[dest]);
        head[dest] = newNode;
    }
}

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

bool Graph::hasEdge(int src, int dest) const
{
    if (src >= numVertices || dest >= numVertices || src < 0 || dest < 0)
        return false;

    Node *current = head[src];
    while (current != nullptr)
    {
        if (current->val == dest)
            return true;
        current = current->next;
    }
    return false;
}

void Graph::addEdge(int src, int dest, int weight)
{
    // Check vertex validity
    if (src >= numVertices || dest >= numVertices || src < 0 || dest < 0)
    {
        std::cout << "Error: Invalid vertex numbers!" << std::endl;
        return;
    }

    // Check for self-loop
    if (src == dest)
    {
        std::cout << "Error: Cannot add self-loop to vertex " << src << std::endl;
        return;
    }

    // Check for duplicate edge
    if (hasEdge(src, dest))
    {
        std::cout << "Error: Edge already exists between " << src << " and " << dest << std::endl;
        return;
    }

    Node *newNode = getAdjacencyList(dest, weight, head[src]);
    head[src] = newNode;

    // undirected graph
    newNode = getAdjacencyList(src, weight, head[dest]);
    head[dest] = newNode;
}

void Graph::removeEdge(int src, int dest)
{
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
            return;
        }
        prev = current;
        current = current->next;
    }
}

void Graph::print_graph(Node *ptr, int i)
{
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