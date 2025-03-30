struct PriorityNode
{
    int data;
    int priority;
    PriorityNode *next;
};

class PriorityQueue
{
private:
    PriorityNode *front;

public:
    PriorityQueue() : front(nullptr) {}
    ~PriorityQueue();
    void enqueue(int value, int priority);
    void dequeue();
    int peek() const;
    bool isEmpty() const;
    void display() const;
};