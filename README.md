# Graph Algorithms Implementation in C++

## Overview
This project implements a graph data structure and several standard graph algorithms in C++. It includes Breadth-First Search (BFS), Depth-First Search (DFS), Dijkstra's shortest path algorithm, and Minimum Spanning Tree algorithms (Prim's and Kruskal's).


## Project Structure
The project is organized into the following components:

### Core Classes
- `Graph`: The main class that represents an undirected weighted graph
- `Algorithms`: A class containing static methods for graph algorithms

### Support Data Structures
- `Queue`: A basic queue implementation for BFS
- `PriorityQueue`: A priority queue implementation for Dijkstra's and Prim's algorithms
- `UnionFind`: A disjoint-set data structure for Kruskal's algorithm

### Files
- `graph.hpp` / `graph.cpp`: Definition and implementation of the Graph class and Algorithms
- `Queue.hpp` / `Queue.cpp`: A simple queue implementation
- `PriorityQueue.hpp` / `PriorityQueue.cpp`: Priority queue implementation
- `UnionFind.hpp` / `UnionFind.cpp`: Union-Find (Disjoint-Set) implementation
- `main.cpp`: Sample program demonstrating the use of the graph algorithms
- `test.cpp`: Comprehensive test suite using doctest

## Features

### Graph Operations
- Create a graph with a specified number of vertices
- Add and remove edges between vertices
- Check if an edge exists between vertices
- Access and print the adjacency list of any vertex

### Graph Algorithms
1. **Breadth-First Search (BFS)**
   - Traverses the graph in breadth-first order
   - Returns a tree representing the BFS traversal

2. **Depth-First Search (DFS)**
   - Traverses the graph in depth-first order
   - Returns a tree representing the DFS traversal

3. **Dijkstra's Algorithm**
   - Finds the shortest paths from a source vertex to all other vertices
   - Returns a tree representing the shortest paths

4. **Prim's Algorithm**
   - Finds a Minimum Spanning Tree (MST) of the graph
   - Uses a priority queue for efficient implementation

5. **Kruskal's Algorithm**
   - Finds a Minimum Spanning Tree (MST) of the graph
   - Uses a Union-Find data structure to detect cycles

## Usage Example
```cpp
#include "graph.hpp"
#include <iostream>

int main() {
    // Create a graph with 5 vertices
    graph::Graph g(5);
    
    // Add edges
    g.addEdge(0, 1, 2);  // Edge from 0 to 1 with weight 2
    g.addEdge(0, 3, 6);  // Edge from 0 to 3 with weight 6
    g.addEdge(1, 2, 3);  // Edge from 1 to 2 with weight 3
    g.addEdge(1, 3, 8);  // Edge from 1 to 3 with weight 8
    g.addEdge(1, 4, 5);  // Edge from 1 to 4 with weight 5
    g.addEdge(2, 4, 7);  // Edge from 2 to 4 with weight 7
    g.addEdge(3, 4, 9);  // Edge from 3 to 4 with weight 9
    
    // Print the graph
    for (int i = 0; i < g.getNumVertices(); i++) {
        g.print_graph(g.getHead(i), i);
    }
    
    // Run BFS starting from vertex 0
    graph::Graph bfsTree = graph::Algorithms::bfs(g, 0);
    
    // Run Prim's algorithm to find MST
    graph::Graph mst = graph::Algorithms::prim(g);
    
    return 0;
}
```