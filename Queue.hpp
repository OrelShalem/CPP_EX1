#include <iostream>
using namespace std;

class Queue
{
private:
    int *arr;
    int front;
    int rear;
    int capacity;
    int count;

public:
    Queue(int size);

    ~Queue();

    void enqueue(int x);

    int dequeue();

    int peek();

    int getSize();

    bool isEmpty();

    bool isFull();
};
