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
    Edge invalidEdges[] = {{0, 5, 1}, {-1, 2, 2}}; // Invalid vertices
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

// בדיקות לאלגוריתם BFS
TEST_CASE("Testing BFS algorithm")
{
    // יצירת גרף מבחן בסיסי
    Edge edges[] = {{0, 1, 1}, {0, 2, 2}, {1, 3, 3}, {2, 3, 4}, {3, 4, 5}};
    Graph g(edges, 5, 5);

    // בדיקת BFS מצומת ההתחלה 0
    Graph bfsResult = Algorithms::bfs(g, 0);

    // בדיקה שהתוצאה היא עץ (ללא מעגלים)
    CHECK(bfsResult.hasEdge(0, 1));
    CHECK(bfsResult.hasEdge(0, 2));

    // בדיקה שרק אחד מהמסלולים נבחר (לא שניהם)
    bool hasEdge1_3 = bfsResult.hasEdge(1, 3);
    bool hasEdge2_3 = bfsResult.hasEdge(2, 3);
    CHECK((hasEdge1_3 || hasEdge2_3));       // לפחות אחד מהם חייב להיות נכון
    CHECK_FALSE((hasEdge1_3 && hasEdge2_3)); // אבל לא שניהם יחד

    CHECK(bfsResult.hasEdge(3, 4));

    // בדיקת מקרי קצה - גרף ריק
    Graph emptyGraph(3);
    Graph emptyBfsResult = Algorithms::bfs(emptyGraph, 0);

    // בגרף ריק, ה-BFS לא יוצר קשתות כלל
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

    // בדיקת מקרה קצה - צומת התחלה שגוי
    CHECK_THROWS_AS(Algorithms::bfs(g, -1), std::invalid_argument); // צומת שלילי
    CHECK_THROWS_AS(Algorithms::bfs(g, 5), std::invalid_argument);  // צומת מחוץ לטווח

    // בדיקת גרף עם רכיב קשירות אחד
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 1);
    disconnectedGraph.addEdge(1, 2, 2);
    // צמתים 3 ו-4 מנותקים

    Graph disconnectedBfs = Algorithms::bfs(disconnectedGraph, 0);

    // בדיקה שרק הצמתים המחוברים נמצאים בתוצאת ה-BFS
    CHECK(disconnectedBfs.hasEdge(0, 1));
    CHECK(disconnectedBfs.hasEdge(1, 2));
    CHECK_FALSE(disconnectedBfs.hasEdge(0, 3));
    CHECK_FALSE(disconnectedBfs.hasEdge(0, 4));
    CHECK_FALSE(disconnectedBfs.hasEdge(1, 3));
    CHECK_FALSE(disconnectedBfs.hasEdge(2, 3));
    CHECK_FALSE(disconnectedBfs.hasEdge(2, 4));

    // בדיקת גרף עם צומת בודד
    Graph singleVertex(1);
    Graph singleBfs = Algorithms::bfs(singleVertex, 0);
    CHECK_EQ(singleBfs.getNumVertices(), 1);

    // בדיקת שצומת ההתחלה בממש קיים ברכיב קשירות אחר
    Graph startFromDisconnected = Algorithms::bfs(disconnectedGraph, 3);
    CHECK_FALSE(startFromDisconnected.hasEdge(3, 0));
    CHECK_FALSE(startFromDisconnected.hasEdge(3, 1));
    CHECK_FALSE(startFromDisconnected.hasEdge(3, 2));
}

// בדיקות לאלגוריתם DFS
TEST_CASE("Testing DFS algorithm")
{
    // יצירת גרף מבחן בסיסי
    Edge edges[] = {{0, 1, 1}, {0, 2, 2}, {1, 3, 3}, {2, 3, 4}, {3, 4, 5}};
    Graph g(edges, 5, 5);

    // בדיקת DFS מצומת ההתחלה 0
    Graph dfsResult = Algorithms::dfs(g, 0);

    // בדיקה שהתוצאה היא עץ (ללא מעגלים)
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
    // כיוון שזה גרף לא מכוון, כל קשת נספרת פעמיים
    CHECK_EQ(edgeCount, (g.getNumVertices() - 1) * 2); // מספר הקשתות בעץ

    // בדיקה שכל הצמתים מחוברים
    bool reachable[5] = {false};

    // סימון הצמתים שניתן להגיע אליהם מצומת 0
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

    // בדיקה שכל הצמתים נגישים
    for (int i = 0; i < 5; i++)
    {
        CHECK(reachable[i]);
    }

    // בדיקת מקרי קצה - גרף ריק
    Graph emptyGraph(3);
    Graph emptyDfsResult = Algorithms::dfs(emptyGraph, 0);

    // בגרף ריק, ה-DFS לא יוצר קשתות כלל
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

    // בדיקת מקרה קצה - צומת התחלה שגוי
    CHECK_THROWS_AS(Algorithms::dfs(g, -1), std::invalid_argument); // צומת שלילי
    CHECK_THROWS_AS(Algorithms::dfs(g, 5), std::invalid_argument);  // צומת מחוץ לטווח

    // בדיקת גרף לא קשיר
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 1);
    disconnectedGraph.addEdge(1, 2, 2);
    disconnectedGraph.addEdge(3, 4, 3); // רכיב קשירות נפרד

    Graph disconnectedDfs = Algorithms::dfs(disconnectedGraph, 0);

    // בדיקה שרק הצמתים בחלק המחובר מופיעים בתוצאה
    CHECK(disconnectedDfs.hasEdge(0, 1));
    CHECK(disconnectedDfs.hasEdge(1, 2));
    CHECK_FALSE(disconnectedDfs.hasEdge(0, 3));
    CHECK_FALSE(disconnectedDfs.hasEdge(0, 4));
    CHECK_FALSE(disconnectedDfs.hasEdge(1, 3));
    CHECK_FALSE(disconnectedDfs.hasEdge(2, 3));
    CHECK_FALSE(disconnectedDfs.hasEdge(3, 4)); // הקשת הזו לא אמורה להיות בתוצאות ה-DFS שמתחיל מ-0

    // בדיקת DFS מהרכיב האחר
    Graph otherComponentDfs = Algorithms::dfs(disconnectedGraph, 3);
    CHECK_FALSE(otherComponentDfs.hasEdge(3, 0));
    CHECK_FALSE(otherComponentDfs.hasEdge(3, 2));
    CHECK(otherComponentDfs.hasEdge(3, 4));
}

// בדיקות לאלגוריתם דייקסטרה
TEST_CASE("Testing Dijkstra algorithm")
{
    // יצירת גרף מבחן עם משקלים
    Edge edges[] = {
        {0, 1, 4},
        {0, 2, 3},
        {1, 2, 2},
        {1, 3, 5},
        {2, 3, 7}};
    Graph g(edges, 5, 4);

    // בדיקת דייקסטרה מצומת ההתחלה 0
    Graph dijkstraResult = Algorithms::dijkstra(g, 0);

    // בדיקה שהתוצאה היא עץ המסלולים הקצרים ביותר
    // מסלול מ-0 ל-1 ישירות
    CHECK(dijkstraResult.hasEdge(0, 1));

    // מסלול מ-0 ל-2 ישירות
    CHECK(dijkstraResult.hasEdge(0, 2));

    // מסלול מ-0 ל-3 דרך 1 או דרך 2
    bool validPathTo3 = (dijkstraResult.hasEdge(1, 3) || dijkstraResult.hasEdge(2, 3));
    CHECK(validPathTo3);

    // בדיקת המשקלים של המסלולים
    if (dijkstraResult.hasEdge(1, 3))
    {
        // אם המסלול ל-3 הוא דרך 1, המשקל של הקשת צריך להיות 5
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
        // אם המסלול ל-3 הוא דרך 2, המשקל של הקשת צריך להיות 7
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

    // בדיקת מקרי קצה - גרף ריק
    Graph emptyGraph(3);
    Graph emptyDijkstraResult = Algorithms::dijkstra(emptyGraph, 0);

    // בגרף ריק, דייקסטרה לא יוצר קשתות כלל
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

    // בדיקת מקרה קצה - צומת התחלה שגוי
    CHECK_THROWS_AS(Algorithms::dijkstra(g, -1), std::invalid_argument); // צומת שלילי
    CHECK_THROWS_AS(Algorithms::dijkstra(g, 5), std::invalid_argument);  // צומת מחוץ לטווח

    // בדיקת גרף לא קשיר
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 10);
    disconnectedGraph.addEdge(1, 2, 20);
    // צמתים 3 ו-4 מנותקים

    Graph disconnectedDijkstra = Algorithms::dijkstra(disconnectedGraph, 0);

    // בדיקה שרק הצמתים המחוברים מופיעים בתוצאה
    CHECK(disconnectedDijkstra.hasEdge(0, 1));
    CHECK(disconnectedDijkstra.hasEdge(1, 2));
    CHECK_FALSE(disconnectedDijkstra.hasEdge(0, 3));
    CHECK_FALSE(disconnectedDijkstra.hasEdge(0, 4));

    // בדיקת גרף עם מסלול רגיל ומסלול "זול" יותר שדורש יותר קשתות
    Graph pathChoiceGraph(4);
    pathChoiceGraph.addEdge(0, 3, 10); // מסלול ישיר יקר
    pathChoiceGraph.addEdge(0, 1, 1);  // מסלול עקיף דרך 1 ו-2
    pathChoiceGraph.addEdge(1, 2, 2);
    pathChoiceGraph.addEdge(2, 3, 3);

    Graph pathChoiceResult = Algorithms::dijkstra(pathChoiceGraph, 0);

    // דייקסטרה צריך לבחור במסלול העקיף כי הוא יותר זול (1+2+3=6 < 10)
    CHECK(pathChoiceResult.hasEdge(0, 1));
    CHECK(pathChoiceResult.hasEdge(1, 2));
    CHECK(pathChoiceResult.hasEdge(2, 3));
    CHECK_FALSE(pathChoiceResult.hasEdge(0, 3)); // הקשת היקרה לא צריכה להיות בעץ
}

// בדיקות לאלגוריתם פרים
TEST_CASE("Testing Prim algorithm")
{
    // יצירת גרף מבחן עם משקלים
    Edge edges[] = {
        {0, 1, 10},
        {0, 2, 6},
        {0, 3, 5},
        {1, 3, 15},
        {2, 3, 4}};
    Graph g(edges, 5, 4);

    // הפעלת אלגוריתם פרים
    Graph primResult = Algorithms::prim(g);

    // בדיקה שהתוצאה היא עץ פורש (n-1 קשתות בגרף עם n צמתים)
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
    // בגרף לא מכוון, כל קשת נספרת פעמיים
    CHECK_EQ(edgeCount, (g.getNumVertices() - 1) * 2);

    // בדיקה שהקשתות שנבחרו הן בעלות המשקל המינימלי
    // העץ הפורש המינימלי צפוי לכלול את הקשתות: (0,3), (2,3), (0,2) עם משקל כולל 5+4+6=15
    CHECK(primResult.hasEdge(0, 3));
    CHECK(primResult.hasEdge(2, 3));

    // אחד מהם אמור להיות נכון
    bool hasEdge0_2 = primResult.hasEdge(0, 2);
    bool hasEdge0_1 = primResult.hasEdge(0, 1);
    CHECK((hasEdge0_2 || hasEdge0_1));

    CHECK_FALSE(primResult.hasEdge(1, 3)); // קשת כבדה שלא אמורה להיות בעץ

    // בדיקת מקרי קצה - גרף ריק
    Graph emptyGraph(3);
    Graph emptyPrimResult = Algorithms::prim(emptyGraph);

    // בגרף ריק, פרים לא יוצר קשתות כלל
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i != j)
            {
                CHECK_FALSE(emptyPrimResult.hasEdge(i, j));
            }
        }
    }

    // בדיקת גרף עם צומת בודד
    Graph singleVertex(1);
    Graph singlePrim = Algorithms::prim(singleVertex);
    CHECK_EQ(singlePrim.getNumVertices(), 1);

    // בדיקת גרף לא קשיר
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 10);
    disconnectedGraph.addEdge(2, 3, 20);
    // שני רכיבים לא מחוברים

    Graph disconnectedPrim = Algorithms::prim(disconnectedGraph);

    // אלגוריתם פרים מתחיל מצומת 0, ולכן רק הרכיב של צומת 0
    // יכלל בעץ הפורש. מכיוון שצמתים 2 ו-3 לא מחוברים לרכיב זה,
    // הם לא יכללו בעץ הפורש לפי המימוש הרגיל של פרים.

    // בודקים שהקשת בין 0 ל-1 נמצאת
    CHECK(disconnectedPrim.hasEdge(0, 1));

    // בודקים שהקשת בין 2 ל-3 לא נמצאת בתוצאה
    CHECK_FALSE(disconnectedPrim.hasEdge(2, 3));

    // לא יהיו קשתות בין רכיבים לא קשירים
    CHECK_FALSE(disconnectedPrim.hasEdge(0, 2)); // אלה שני רכיבים נפרדים
    CHECK_FALSE(disconnectedPrim.hasEdge(1, 3));
}

// בדיקות לאלגוריתם קרוסקל
TEST_CASE("Testing Kruskal algorithm")
{
    // יצירת גרף מבחן עם משקלים
    Edge edges[] = {
        {0, 1, 10},
        {0, 2, 6},
        {0, 3, 5},
        {1, 3, 15},
        {2, 3, 4}};
    Graph g(edges, 5, 4);

    // הפעלת אלגוריתם קרוסקל
    Graph kruskalResult = Algorithms::kruskal(g);

    // בדיקה שהתוצאה היא עץ פורש (n-1 קשתות בגרף עם n צמתים)
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
    // בגרף לא מכוון, כל קשת נספרת פעמיים
    CHECK_EQ(edgeCount, (g.getNumVertices() - 1) * 2);

    // בדיקה שהקשתות שנבחרו הן בעלות המשקל המינימלי
    // העץ הפורש המינימלי צפוי לכלול את הקשתות: (2,3), (0,3), (0,2) עם משקל כולל 4+5+6=15
    CHECK(kruskalResult.hasEdge(2, 3)); // קשת עם המשקל הקטן ביותר
    CHECK(kruskalResult.hasEdge(0, 3)); // קשת עם המשקל השני הקטן ביותר

    // אחד מהם אמור להיות נכון
    bool hasEdge0_2 = kruskalResult.hasEdge(0, 2);
    bool hasEdge0_1 = kruskalResult.hasEdge(0, 1);
    CHECK((hasEdge0_2 || hasEdge0_1));

    CHECK_FALSE(kruskalResult.hasEdge(1, 3)); // קשת כבדה שלא אמורה להיות בעץ

    // בדיקת מקרי קצה - גרף ריק
    Graph emptyGraph(3);
    Graph emptyKruskalResult = Algorithms::kruskal(emptyGraph);

    // בגרף ריק, קרוסקל לא יוצר קשתות כלל
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i != j)
            {
                CHECK_FALSE(emptyKruskalResult.hasEdge(i, j));
            }
        }
    }

    // בדיקת גרף עם צומת בודד
    Graph singleVertex(1);
    Graph singleKruskal = Algorithms::kruskal(singleVertex);
    CHECK_EQ(singleKruskal.getNumVertices(), 1);

    // בדיקת גרף לא קשיר
    Graph disconnectedGraph(5);
    disconnectedGraph.addEdge(0, 1, 10);
    disconnectedGraph.addEdge(2, 3, 20);
    disconnectedGraph.addEdge(3, 4, 30);
    // שני רכיבים לא מחוברים

    Graph disconnectedKruskal = Algorithms::kruskal(disconnectedGraph);

    // קרוסקל צריך ליצור יער (מספר עצים פורשים, אחד לכל רכיב קשירות)
    CHECK(disconnectedKruskal.hasEdge(0, 1));
    CHECK(disconnectedKruskal.hasEdge(2, 3));
    CHECK(disconnectedKruskal.hasEdge(3, 4));
    CHECK_FALSE(disconnectedKruskal.hasEdge(0, 2)); // אלה שני רכיבים נפרדים
    CHECK_FALSE(disconnectedKruskal.hasEdge(1, 4));

    // בדיקת גרף עם מספר קשתות באותו משקל
    Graph equalWeightGraph(4);
    equalWeightGraph.addEdge(0, 1, 5);
    equalWeightGraph.addEdge(0, 2, 5);
    equalWeightGraph.addEdge(0, 3, 5);
    equalWeightGraph.addEdge(1, 2, 5);
    equalWeightGraph.addEdge(2, 3, 5);

    Graph equalWeightResult = Algorithms::kruskal(equalWeightGraph);

    // בדיקה שהתוצאה היא עץ פורש (n-1 קשתות)
    edgeCount = 0;
    for (int i = 0; i < equalWeightGraph.getNumVertices(); i++)
    {
        Node *head = equalWeightResult.getHead(i);
        while (head != nullptr)
        {
            edgeCount++;
            head = head->next;
        }
    }
    CHECK_EQ(edgeCount, (equalWeightGraph.getNumVertices() - 1) * 2);
}

// בדיקות למבנה נתונים Queue
TEST_CASE("Testing Queue")
{
    // בדיקת יצירת תור
    Queue q(5);
    CHECK(q.isEmpty());
    CHECK_FALSE(q.isFull());
    CHECK_EQ(q.getSize(), 0);

    // בדיקת הכנסה לתור והוצאה מתור
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

    // בדיקת מקרי קצה
    CHECK_THROWS_AS(q.peek(), std::underflow_error);    // תור ריק
    CHECK_THROWS_AS(q.dequeue(), std::underflow_error); // תור ריק

    // מילוי התור
    q.enqueue(1);
    q.enqueue(2);
    q.enqueue(3);
    q.enqueue(4);
    q.enqueue(5);
    CHECK(q.isFull());

    // ניסיון להכניס לתור מלא
    CHECK_THROWS_AS(q.enqueue(6), std::overflow_error);

    // בדיקת התנהגות מעגלית
    CHECK_EQ(q.dequeue(), 1);
    CHECK_EQ(q.dequeue(), 2);

    // עכשיו התור לא מלא
    CHECK_FALSE(q.isFull());

    // ניתן להכניס שוב
    q.enqueue(6);
    q.enqueue(7);
    CHECK(q.isFull());

    // בדיקה שהתור עובד בצורה מעגלית
    CHECK_EQ(q.dequeue(), 3);
    CHECK_EQ(q.dequeue(), 4);
    CHECK_EQ(q.dequeue(), 5);
    CHECK_EQ(q.dequeue(), 6);
    CHECK_EQ(q.dequeue(), 7);
    CHECK(q.isEmpty());
}

// בדיקות למבנה נתונים PriorityQueue
TEST_CASE("Testing PriorityQueue")
{
    // בדיקת יצירת תור עדיפויות
    PriorityQueue pq;
    CHECK(pq.isEmpty());

    // בדיקת הכנסה לתור והוצאה לפי עדיפות
    pq.enqueue(10, 2); // ערך 10 עם עדיפות 2
    pq.enqueue(20, 1); // ערך 20 עם עדיפות 1 (גבוהה יותר)
    pq.enqueue(30, 3); // ערך 30 עם עדיפות 3 (נמוכה יותר)

    // בדיקה שפעולת peek מחזירה את הערך עם העדיפות הגבוהה ביותר
    CHECK_EQ(pq.peek(), 20);

    // הוצאה מהתור
    pq.dequeue();
    CHECK_EQ(pq.peek(), 10);

    pq.dequeue();
    CHECK_EQ(pq.peek(), 30);

    pq.dequeue();
    CHECK(pq.isEmpty());

    // בדיקת מקרי קצה
    CHECK_THROWS_AS(pq.peek(), std::underflow_error);    // תור ריק
    CHECK_THROWS_AS(pq.dequeue(), std::underflow_error); // תור ריק

    // בדיקת תור עם עדיפויות זהות
    pq.enqueue(100, 5);
    pq.enqueue(200, 5);
    pq.enqueue(300, 5);

    // בדיקה שהערכים יוצאים לפי סדר הכנסתם כאשר העדיפויות זהות
    CHECK_EQ(pq.peek(), 100);
    pq.dequeue();
    CHECK_EQ(pq.peek(), 200);
    pq.dequeue();
    CHECK_EQ(pq.peek(), 300);
    pq.dequeue();
    CHECK(pq.isEmpty());

    // בדיקת עדיפויות מעורבות
    pq.enqueue(1, 10);
    pq.enqueue(2, 5);
    pq.enqueue(3, 15);
    pq.enqueue(4, 7);
    pq.enqueue(5, 3);

    // בדיקה שהערכים יוצאים בסדר עדיפויות עולה
    CHECK_EQ(pq.peek(), 5); // עדיפות 3
    pq.dequeue();
    CHECK_EQ(pq.peek(), 2); // עדיפות 5
    pq.dequeue();
    CHECK_EQ(pq.peek(), 4); // עדיפות 7
    pq.dequeue();
    CHECK_EQ(pq.peek(), 1); // עדיפות 10
    pq.dequeue();
    CHECK_EQ(pq.peek(), 3); // עדיפות 15
    pq.dequeue();
    CHECK(pq.isEmpty());
}

// בדיקות למבנה נתונים UnionFind
TEST_CASE("Testing UnionFind")
{
    // בדיקת יצירת מבנה Union-Find
    UnionFind uf(5);

    // בבדיקת במצב התחלתי כל צומת הוא בקבוצה משלו
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

    // בדיקה שיש בדיוק 5 קבוצות נפרדות
    CHECK_EQ(uf.getSetCount(), 5);

    // בדיקת איחוד צמתים
    uf.unionSet(0, 1);
    CHECK(uf.connected(0, 1));
    CHECK_EQ(uf.getSetCount(), 4);

    uf.unionSet(2, 3);
    CHECK(uf.connected(2, 3));
    CHECK_EQ(uf.getSetCount(), 3);

    // בדיקה שהצמתים עדיין לא מחוברים
    CHECK_FALSE(uf.connected(0, 2));
    CHECK_FALSE(uf.connected(1, 3));

    // בדיקה שאיחוד נוסף עובד בצורה טרנזיטיבית
    uf.unionSet(0, 2);
    CHECK(uf.connected(0, 2));
    CHECK(uf.connected(0, 3)); // 0 מחובר ל-2 ו-2 מחובר ל-3
    CHECK(uf.connected(1, 3)); // 1 מחובר ל-0, ו-0 מחובר ל-3
    CHECK_EQ(uf.getSetCount(), 2);

    // בדיקה שצומת 4 עדיין מנותק
    CHECK_FALSE(uf.connected(0, 4));
    CHECK_FALSE(uf.connected(1, 4));
    CHECK_FALSE(uf.connected(2, 4));
    CHECK_FALSE(uf.connected(3, 4));

    // איחוד אחרון
    uf.unionSet(3, 4);
    CHECK(uf.connected(0, 4)); // כל הצמתים מחוברים כעת
    CHECK_EQ(uf.getSetCount(), 1);

    // בדיקת מקרי קצה - איחוד קבוצה עם עצמה
    int setCountBefore = uf.getSetCount();
    uf.unionSet(0, 0);
    CHECK_EQ(uf.getSetCount(), setCountBefore); // לא השתנה

    // ניסיון לאחד צמתים שכבר באותה קבוצה
    setCountBefore = uf.getSetCount();
    uf.unionSet(1, 4);
    CHECK_EQ(uf.getSetCount(), setCountBefore); // לא השתנה
}