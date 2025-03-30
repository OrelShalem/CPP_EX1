#include "graph.hpp"
#include <iostream>
#include <climits>
#include "Queue.hpp"
#include "PriorityQueue.hpp"
#include "UnionFind.hpp"

namespace graph
{

    // constructor
    Graph::Graph(Edge edges[], int n, int numVertices) : numVertices(numVertices)
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
            // 0 -> 1
            Node *newNode = getAdjacencyList(dest, weight, head[src]);
            head[src] = newNode;

            // undirected graph
            // 1 -> 0
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
    }

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
    }

    // get number of vertices in the graph (constant)
    int Graph::getNumVertices() const
    {
        return numVertices;
    }

    Node *Graph::getHead(int src) const
    {
        if (src >= getNumVertices() || src < 0)
            throw std::invalid_argument("Invalid vertex number!");
        return head[src];
    }

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

    // מימוש אלגוריתם BFS
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

        // יצירת תור לאלגוריתם BFS
        Queue q(numVertices);

        // מסמנים את צומת ההתחלה כמבוקר ומכניסים לתור
        visited[start] = true;
        q.enqueue(start);

        // מבצעים BFS
        while (!q.isEmpty())
        {
            // מוציאים צומת מהתור
            int current = q.dequeue();

            // עוברים על כל השכנים של הצומת הנוכחי
            Node *neighbor = graph.getHead(current);
            while (neighbor != nullptr)
            {
                int adjVertex = neighbor->val;
                int weight = neighbor->cost;

                // אם לא ביקרנו בשכן זה עדיין
                if (!visited[adjVertex])
                {
                    // מוסיפים קשת לגרף התוצאה
                    result.addEdge(current, adjVertex, weight);

                    // מסמנים את השכן כמבוקר ומכניסים לתור
                    visited[adjVertex] = true;
                    q.enqueue(adjVertex);
                }

                // ממשיכים לשכן הבא
                neighbor = neighbor->next;
            }
        }

        // שחרור זיכרון
        delete[] visited;

        return result;
    }

    // פונקציית עזר לאלגוריתם DFS
    void dfsUtil(const Graph &graph, int vertex, bool *visited, Graph &result, int parent)
    {
        // מסמנים את הצומת הנוכחי כמבוקר
        visited[vertex] = true;

        // עוברים על כל השכנים של הצומת הנוכחי
        Node *neighbor = graph.getHead(vertex);
        while (neighbor != nullptr)
        {
            int adjVertex = neighbor->val;
            int weight = neighbor->cost;

            // אם לא ביקרנו בשכן זה עדיין
            if (!visited[adjVertex])
            {
                // מוסיפים קשת לגרף התוצאה
                result.addEdge(vertex, adjVertex, weight);

                // קריאה רקורסיבית לשכן
                dfsUtil(graph, adjVertex, visited, result, vertex);
            }
            // ניתן להוסיף בדיקה כאן אם רוצים לטפל בקשתות אחורה
            // בדיקה זו יכולה להשתמש בפרמטר parent כדי להימנע מספירת קשת הדדית כמו קשת אחורה
            else if (adjVertex != parent)
            {
                // קשת אחורה נמצאה (לא משנים את גרף התוצאה עבור קשתות אחורה)
                // std::cout << "Back edge from " << vertex << " to " << adjVertex << std::endl;
            }

            // ממשיכים לשכן הבא
            neighbor = neighbor->next;
        }
    }

    // מימוש אלגוריתם DFS
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

    // מימוש אלגוריתם דייקסטרה
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

            // מסמנים את הצומת כמסתיים
            finalized[u] = true;

            // עוברים על כל השכנים של הצומת הנוכחי
            Node *neighbor = graph.getHead(u);
            while (neighbor != nullptr)
            {
                int v = neighbor->val;
                int weight = neighbor->cost;

                // אם מצאנו מסלול קצר יותר
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

        // בניית גרף התוצאה מעץ המסלולים הקצרים ביותר
        for (int i = 0; i < numVertices; i++)
        {
            if (i != start && parent[i] != -1 && distance[i] != INT_MAX)
            {
                // מוסיפים קשת מההורה לצומת הנוכחי
                // המשקל של הקשת הוא ההפרש במרחקים
                int edgeWeight = distance[i] - distance[parent[i]];
                result.addEdge(parent[i], i, edgeWeight);
            }
        }

        // שחרור זיכרון
        delete[] distance;
        delete[] finalized;
        delete[] parent;

        return result;
    }

    // מימוש אלגוריתם פרים
    Graph Algorithms::prim(const Graph &graph)
    {
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

    // מימוש אלגוריתם קרוסקל
    Graph Algorithms::kruskal(const Graph &graph)
    {
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
}