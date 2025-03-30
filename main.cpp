#include <iostream>
#include <iomanip>
#include "graph.hpp"

using namespace graph;

// פונקציה להדפסת כותרת מודגשת
void printHeader(const std::string &title)
{
    std::cout << "\n\n"
              << std::string(50, '=') << std::endl;
    std::cout << std::setw(25 + title.length() / 2) << title << std::endl;
    std::cout << std::string(50, '=') << std::endl;
}

// יצירת גרף לדוגמה באופן אוטומטי - ללא צורך בקלט מהמשתמש
Graph createExampleGraph()
{
    // יצירת גרף עם 6 צמתים
    int numVertices = 6;
    int numEdges = 9;

    Edge edges[9] = {
        {0, 1, 4}, // קשת מ-0 ל-1 במשקל 4
        {0, 2, 3}, // קשת מ-0 ל-2 במשקל 3
        {1, 2, 2}, // קשת מ-1 ל-2 במשקל 2
        {1, 3, 5}, // קשת מ-1 ל-3 במשקל 5
        {2, 3, 7}, // קשת מ-2 ל-3 במשקל 7
        {2, 4, 8}, // קשת מ-2 ל-4 במשקל 8
        {3, 4, 1}, // קשת מ-3 ל-4 במשקל 1
        {3, 5, 6}, // קשת מ-3 ל-5 במשקל 6
        {4, 5, 9}  // קשת מ-4 ל-5 במשקל 9
    };

    return Graph(edges, numEdges, numVertices);
}

// פונקציה המאפשרת למשתמש ליצור גרף באמצעות קלט
Graph createUserGraph()
{
    int numVertices;

    std::cout << "Enter number of vertices: ";
    std::cin >> numVertices;

    std::cout << "Enter number of edges: ";
    int numEdges;
    std::cin >> numEdges;

    Edge *edges = new Edge[numEdges];
    int validEdges = 0; // מעקב אחר מספר הקשתות התקינות

    std::cout << "\nFor each edge enter: source vertex, destination vertex, weight\n";
    for (int i = 0; i < numEdges; i++)
    {
        std::cout << "Edge " << i + 1 << ": ";
        std::cin >> edges[i].src >> edges[i].dest >> edges[i].weight;

        // בדיקת תקינות - בדיקת מספרי הצמתים
        if (edges[i].src >= numVertices || edges[i].dest >= numVertices ||
            edges[i].src < 0 || edges[i].dest < 0)
        {
            std::cout << "Invalid vertex numbers! Vertices should be between 0 and "
                      << numVertices - 1 << std::endl;
            i--; // נסה שוב עם אותו מיקום במערך
            continue;
        }

        // בדיקת לולאות עצמיות
        if (edges[i].src == edges[i].dest)
        {
            std::cout << "Self-loops are not allowed in a simple graph! Please try again." << std::endl;
            i--; // נסה שוב עם אותו מיקום במערך
            continue;
        }

        // בדיקת קשתות כפולות
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
            i--; // נסה שוב עם אותו מיקום במערך
            continue;
        }

        // אם הגענו לכאן, הקשת תקינה
        validEdges++;
    }

    // יצירת הגרף
    Graph graph(edges, validEdges, numVertices);
    delete[] edges;

    return graph;
}

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

    Graph graph = (choice == 1) ? createExampleGraph() : createUserGraph();

    // הדפסת הגרף המקורי
    printHeader("ORIGINAL GRAPH");
    std::cout << "Graph structure with " << graph.getNumVertices() << " vertices:" << std::endl;
    for (int i = 0; i < graph.getNumVertices(); i++)
    {
        graph.print_graph(graph.getHead(i), i);
    }

    // בחירת צומת התחלה לאלגוריתמים
    int startVertex;
    std::cout << "\nChoose a start vertex for BFS, DFS and Dijkstra (0-" << graph.getNumVertices() - 1 << "): ";
    std::cin >> startVertex;

    // בדיקת תקינות צומת ההתחלה
    if (startVertex < 0 || startVertex >= graph.getNumVertices())
    {
        std::cout << "Invalid vertex! Using vertex 0 instead." << std::endl;
        startVertex = 0;
    }

    try
    {
        // אלגוריתם BFS
        printHeader("BFS ALGORITHM");
        std::cout << "Running BFS starting from vertex " << startVertex << ":" << std::endl;
        Graph bfsResult = Algorithms::bfs(graph, startVertex);
        std::cout << "BFS result (tree structure):" << std::endl;
        for (int i = 0; i < bfsResult.getNumVertices(); i++)
        {
            bfsResult.print_graph(bfsResult.getHead(i), i);
        }

        // אלגוריתם DFS
        printHeader("DFS ALGORITHM");
        std::cout << "Running DFS starting from vertex " << startVertex << ":" << std::endl;
        Graph dfsResult = Algorithms::dfs(graph, startVertex);
        std::cout << "DFS result (tree structure):" << std::endl;
        for (int i = 0; i < dfsResult.getNumVertices(); i++)
        {
            dfsResult.print_graph(dfsResult.getHead(i), i);
        }

        // אלגוריתם דייקסטרה
        printHeader("DIJKSTRA ALGORITHM");
        std::cout << "Running Dijkstra's algorithm from vertex " << startVertex << ":" << std::endl;
        Graph dijkstraResult = Algorithms::dijkstra(graph, startVertex);
        std::cout << "Dijkstra result (shortest paths tree):" << std::endl;
        for (int i = 0; i < dijkstraResult.getNumVertices(); i++)
        {
            dijkstraResult.print_graph(dijkstraResult.getHead(i), i);
        }

        // אלגוריתם פרים
        printHeader("PRIM'S ALGORITHM");
        std::cout << "Running Prim's algorithm for Minimum Spanning Tree:" << std::endl;
        Graph primResult = Algorithms::prim(graph);
        std::cout << "Prim result (MST):" << std::endl;
        for (int i = 0; i < primResult.getNumVertices(); i++)
        {
            primResult.print_graph(primResult.getHead(i), i);
        }

        // אלגוריתם קרוסקל
        printHeader("KRUSKAL'S ALGORITHM");
        std::cout << "Running Kruskal's algorithm for Minimum Spanning Tree:" << std::endl;
        Graph kruskalResult = Algorithms::kruskal(graph);
        std::cout << "Kruskal result (MST):" << std::endl;
        for (int i = 0; i < kruskalResult.getNumVertices(); i++)
        {
            kruskalResult.print_graph(kruskalResult.getHead(i), i);
        }

        // בדיקת הוספה והסרה של קשתות
        printHeader("EDGE OPERATIONS");
        std::cout << "Demonstrating edge operations:" << std::endl;

        // העתקת הגרף המקורי
        Graph testGraph = createExampleGraph();

        // הדפסת הגרף המקורי
        std::cout << "Original graph:" << std::endl;
        for (int i = 0; i < testGraph.getNumVertices(); i++)
        {
            testGraph.print_graph(testGraph.getHead(i), i);
        }

        // הוספת קשת חדשה
        std::cout << "\nAdding edge from 0 to 5 with weight 10:" << std::endl;
        testGraph.addEdge(0, 5, 10);
        std::cout << "After adding edge:" << std::endl;
        for (int i = 0; i < testGraph.getNumVertices(); i++)
        {
            testGraph.print_graph(testGraph.getHead(i), i);
        }

        // הסרת קשת
        std::cout << "\nRemoving edge between 0 and 1:" << std::endl;
        testGraph.removeEdge(0, 1);
        std::cout << "After removing edge:" << std::endl;
        for (int i = 0; i < testGraph.getNumVertices(); i++)
        {
            testGraph.print_graph(testGraph.getHead(i), i);
        }

        // בדיקת קיום קשת
        std::cout << "\nChecking if edges exist:" << std::endl;
        std::cout << "Edge 0-2 exists: " << (testGraph.hasEdge(0, 2) ? "Yes" : "No") << std::endl;
        std::cout << "Edge 0-1 exists: " << (testGraph.hasEdge(0, 1) ? "Yes" : "No") << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    std::cout << "\nDemonstration completed successfully!" << std::endl;
    return 0;
}