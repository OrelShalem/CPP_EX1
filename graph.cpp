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
                // שימוש ב-addEdge במקום להוסיף ישירות לרשימה
                addEdge(src, dest, weight);
            }
            catch (const std::invalid_argument &e)
            {
                // אם addEdge זורק חריגה, נדפיס אזהרה ונמשיך
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

        // בדיקת תקינות צומת ההתחלה
        if (start < 0 || start >= numVertices)
        {
            throw std::invalid_argument("Invalid start vertex!");
        }

        // יצירת גרף תוצאה ריק
        Graph result(numVertices);

        // מערך לסימון צמתים שכבר ביקרנו בהם
        bool *visited = new bool[numVertices]();

        // תור לאחסון הצמתים שנגלו אך טרם נבדקו
        Queue q(numVertices);

        // סימון צומת ההתחלה כמבוקר
        visited[start] = true;

        // הכנסת צומת ההתחלה לתור
        q.enqueue(start);

        // כל עוד התור לא ריק
        while (!q.isEmpty())
        {
            // הוצאת צומת מהתור
            int u = q.dequeue();

            // קבלת רשימת השכנים של הצומת
            Node *neighbor = graph.getHead(u);

            // עבור כל שכן
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // אם השכן טרם בוקר
                if (!visited[v])
                {
                    // סימון השכן כמבוקר
                    visited[v] = true;

                    // הוספת קשת לגרף התוצאה
                    result.addEdge(u, v, weight);

                    // הכנסת השכן לתור
                    q.enqueue(v);
                }

                // מעבר לשכן הבא
                neighbor = neighbor->next;
            }
        }

        // שחרור זיכרון
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
        // סימון הצומת הנוכחי כמבוקר
        visited[u] = true;

        // קבלת רשימת השכנים של הצומת
        Node *neighbor = graph.getHead(u);

        // עבור כל שכן
        while (neighbor != nullptr)
        {
            int v = neighbor->val;
            int weight = neighbor->cost;

            // אם השכן טרם בוקר
            if (!visited[v])
            {
                // הוספת קשת לגרף התוצאה
                result.addEdge(u, v, weight);

                // קריאה רקורסיבית לשכן - מעביר את הצומת הנוכחי כהורה
                dfsUtil(graph, v, visited, result, u);
            }
            // מונע הוספת קשת חזרה להורה (אם זהו גרף לא מכוון)
            // אם השכן הוא ההורה של הצומת הנוכחי, נדלג עליו
            else if (v != parent)
            {
                // כאן ניתן להוסיף לוגיקה נוספת אם נרצה לטפל בקשתות אחורה
                // לדוגמה: זיהוי מעגלים
            }

            // מעבר לשכן הבא
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

        // בדיקת תקינות צומת ההתחלה
        if (start < 0 || start >= numVertices)
        {
            throw std::invalid_argument("Invalid start vertex!");
        }

        // יצירת גרף תוצאה ריק
        Graph result(numVertices);

        // מערך לסימון צמתים שכבר ביקרנו בהם
        bool *visited = new bool[numVertices]();

        // מבצעים DFS מצומת ההתחלה
        dfsUtil(graph, start, visited, result, -1);

        // שחרור זיכרון
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

        // בדיקת תקינות צומת ההתחלה
        if (start < 0 || start >= numVertices)
        {
            throw std::invalid_argument("Invalid start vertex!");
        }

        // יצירת גרף תוצאה ריק
        Graph result(numVertices);

        // מערך המרחקים הקצרים ביותר מצומת ההתחלה
        int *distance = new int[numVertices];

        // מערך לסימון צמתים שכבר הסתיימו
        bool *finalized = new bool[numVertices]();

        // מערך שומר את ההורה של כל צומת במסלול הקצר ביותר
        int *parent = new int[numVertices];

        // אתחול מערכים
        for (int i = 0; i < numVertices; i++)
        {
            distance[i] = INT_MAX; // אינסוף
            parent[i] = -1;        // אין הורה
        }

        // המרחק מצומת ההתחלה לעצמו הוא 0
        distance[start] = 0;

        // תור עדיפויות למציאת הצומת הבא עם המרחק המינימלי
        PriorityQueue pq;

        // מכניסים את צומת ההתחלה לתור העדיפויות
        pq.enqueue(start, 0);

        // כל עוד התור לא ריק
        while (!pq.isEmpty())
        {
            // מוציאים את הצומת עם המרחק המינימלי
            int u = pq.peek();
            pq.dequeue();

            // אם כבר סיימנו עם צומת זה, נמשיך
            if (finalized[u])
                continue;

            // מסמנים שסיימנו עם צומת זה
            finalized[u] = true;

            // אם יש להורה, מוסיפים את הקשת לגרף התוצאה
            if (parent[u] != -1)
            {
                // המשקל הוא מרחק הצומת ולא משקל הקשת המקורי!
                result.addEdge(parent[u], u, distance[u] - distance[parent[u]]);
            }

            // עוברים על כל השכנים של הצומת הנוכחי
            Node *neighbor = graph.getHead(u);
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // אם השכן טרם סיומי ויש מסלול קצר יותר דרך צומת u
                if (!finalized[v] && distance[u] != INT_MAX &&
                    distance[u] + weight < distance[v])
                {
                    // עדכון המרחק
                    distance[v] = distance[u] + weight;

                    // עדכון ההורה
                    parent[v] = u;

                    // הוספה לתור העדיפויות
                    pq.enqueue(v, distance[v]);
                }

                // ממשיכים לשכן הבא
                neighbor = neighbor->next;
            }
        }

        // שחרור זיכרון
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

        // יצירת גרף תוצאה ריק
        Graph result(numVertices);

        // מערך לסימון צמתים שכבר נכללו בעץ
        bool *inMST = new bool[numVertices]();

        // מערך המשקלים המינימליים לחיבור כל צומת לעץ
        int *key = new int[numVertices];

        // מערך ההורים של כל צומת בעץ הפורש
        int *parent = new int[numVertices];

        // אתחול מערכים
        for (int i = 0; i < numVertices; i++)
        {
            key[i] = INT_MAX; // אינסוף
            parent[i] = -1;   // אין הורה
        }

        // נתחיל מהצומת הראשון (0)
        key[0] = 0; // המשקל לצומת ההתחלה הוא 0

        // תור עדיפויות למציאת הצומת הבא עם המשקל המינימלי
        PriorityQueue pq;

        // מכניסים את צומת ההתחלה לתור העדיפויות
        pq.enqueue(0, 0);

        // כל עוד לא הוספנו את כל הצמתים לעץ
        while (!pq.isEmpty())
        {
            // מוציאים את הצומת עם המשקל המינימלי
            int u = pq.peek();
            pq.dequeue();

            // אם כבר הוספנו צומת זה לעץ, נמשיך
            if (inMST[u])
                continue;

            // מוסיפים את הצומת לעץ
            inMST[u] = true;

            // אם יש להורה, מוסיפים את הקשת לגרף התוצאה
            if (parent[u] != -1)
            {
                result.addEdge(u, parent[u], key[u]);
            }

            // עוברים על כל השכנים של הצומת הנוכחי
            Node *neighbor = graph.getHead(u);
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // אם השכן עדיין לא בעץ ויש קשת קלה יותר
                if (!inMST[v] && weight < key[v])
                {
                    // עדכון המשקל
                    key[v] = weight;

                    // עדכון ההורה
                    parent[v] = u;

                    // הוספה לתור העדיפויות
                    pq.enqueue(v, key[v]);
                }

                // ממשיכים לשכן הבא
                neighbor = neighbor->next;
            }
        }

        // שחרור זיכרון
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

        // יצירת גרף תוצאה ריק
        Graph result(numVertices);

        // איסוף כל הקשתות מהגרף
        // במקרה הגרוע יש n(n-1)/2 קשתות בגרף לא מכוון
        int maxEdges = numVertices * (numVertices - 1) / 2;
        Edge *edges = new Edge[maxEdges];
        int edgeCount = 0;

        // עוברים על כל הצמתים והקשתות
        for (int i = 0; i < numVertices; i++)
        {
            Node *neighbor = graph.getHead(i);
            while (neighbor != nullptr)
            {
                // למנוע כפילויות (כי זה גרף לא מכוון)
                // נוסיף את הקשת רק אם i < שכן
                if (i < neighbor->val)
                {
                    edges[edgeCount].src = i;
                    edges[edgeCount].dest = neighbor->val;
                    edges[edgeCount].weight = neighbor->cost;
                    edgeCount++;
                }

                // ממשיכים לשכן הבא
                neighbor = neighbor->next;
            }
        }

        // מיון הקשתות לפי משקל בסדר עולה
        for (int i = 0; i < edgeCount - 1; i++)
        {
            for (int j = 0; j < edgeCount - i - 1; j++)
            {
                if (edges[j].weight > edges[j + 1].weight)
                {
                    // החלפה
                    Edge temp = edges[j];
                    edges[j] = edges[j + 1];
                    edges[j + 1] = temp;
                }
            }
        }

        // יצירת מבנה Union-Find
        UnionFind uf(numVertices);

        // עוברים על הקשתות בסדר עולה של משקל
        for (int i = 0; i < edgeCount; i++)
        {
            int src = edges[i].src;
            int dest = edges[i].dest;

            // בודקים אם הקשת יוצרת מעגל
            if (!uf.connected(src, dest))
            {
                // אם לא, מוסיפים אותה לעץ הפורש המינימלי
                result.addEdge(src, dest, edges[i].weight);

                // מאחדים את הקבוצות
                uf.unionSet(src, dest);
            }
        }

        // שחרור זיכרון
        delete[] edges;

        return result;
    }

    // bool Algorithms::isConnected(const Graph &graph)
    // {
    //     int numVertices = graph.getNumVertices();

    //     // יצירת גרף תוצאה ריק
    //     Graph result(numVertices);

    //     // מערך לסימון צמתים שכבר ביקרנו בהם
    //     bool *visited = new bool[numVertices]();

    //     // מבצעים DFS מצומת ההתחלה
    //     dfsUtil(graph, 0, visited, result, -1);

    //     // בדיקה אם כל הצמתים נכללו בעץ
    //     for (int i = 0; i < numVertices; i++)
    //     {
    //         if (!visited[i])
    //         {
    //             std::cout << "Graph is not connected!" << std::endl;
    //             // שחרור זיכרון
    //             delete[] visited;
    //             return false;
    //         }
    //     }

    //     // שחרור זיכרון
    //     delete[] visited;

    //     std::cout << "Graph is connected!" << std::endl;
    //     return true;
    // }
}