#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "graph.hpp"

using namespace graph;

TEST_CASE("Graph creation test")
{
    Edge edges[] = {{0, 1, 1}, {1, 2, 2}};
    Graph g(edges, 2, 3);

    CHECK(g.hasEdge(0, 1));
    CHECK(g.hasEdge(1, 2));
    CHECK_FALSE(g.hasEdge(0, 2));
}

TEST_CASE("Edge addition test")
{
    Edge edges[] = {{0, 1, 1}};
    Graph g(edges, 1, 3);

    g.addEdge(1, 2, 2);

    CHECK(g.hasEdge(1, 2));
}

TEST_CASE("Edge removal test")
{
    Edge edges[] = {{0, 1, 1}, {1, 2, 2}};
    Graph g(edges, 2, 3);

    g.removeEdge(0, 1);

    CHECK_FALSE(g.hasEdge(0, 1));
    CHECK(g.hasEdge(1, 2));
}

TEST_CASE("Simple graph test - no self loops or parallel edges")
{
    Edge edges[] = {{0, 1, 1}};
    Graph g(edges, 1, 3);

    // Attempt to create self loop
    g.addEdge(0, 0, 5);
    CHECK_FALSE(g.hasEdge(0, 0));

    // Attempt to create parallel edge
    g.addEdge(0, 1, 10);

    // Check no additional edge was created and weight remained unchanged
    int edgeCount = 0;
    int weight = 0;

    Node *current = g.head[0];
    while (current)
    {
        if (current->val == 1)
        {
            edgeCount++;
            weight = current->cost;
        }
        current = current->next;
    }

    CHECK(edgeCount == 1); // Should only have one edge
    CHECK(weight == 1);    // Weight should remain 1, not 10
}