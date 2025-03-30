#include "Queue.hpp"
#include <iostream>

Queue::Queue(int size) : arr(new int[size]), front(0), rear(-1), capacity(size), count(0) {}

Queue::~Queue()
{
    delete[] arr;
}

void Queue::enqueue(int item)
{
    if (isFull())
    {
        throw std::overflow_error("Queue is full");
    }
    rear = (rear + 1) % capacity;
    arr[rear] = item;
    count++;
}

int Queue::dequeue()
{
    if (isEmpty())
    {
        throw std::underflow_error("Queue is empty");
    }
    int item = arr[front];
    front = (front + 1) % capacity;
    count--;
    return item;
}

int Queue::peek()
{
    if (isEmpty())
    {
        throw std::underflow_error("Queue is empty");
    }
    return arr[front];
}

bool Queue::isEmpty()
{
    return count == 0;
}

bool Queue::isFull()
{
    return count == capacity;
}

int Queue::getSize()
{
    return count;
}