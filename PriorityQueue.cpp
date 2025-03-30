#include "PriorityQueue.hpp"
#include <stdexcept>
#include <iostream>

PriorityQueue::~PriorityQueue()
{
    while (!isEmpty())
    {
        dequeue();
    }
}

void PriorityQueue::enqueue(int value, int priority)
{
    PriorityNode *newNode = new PriorityNode;
    newNode->data = value;
    newNode->priority = priority;
    newNode->next = nullptr;

    if (isEmpty() || priority < front->priority)
    {
        newNode->next = front;
        front = newNode;
    }
    else
    {
        PriorityNode *current = front;
        while (current->next != nullptr && current->next->priority <= priority)
        {
            current = current->next;
        }
        newNode->next = current->next;
        current->next = newNode;
    }
}

void PriorityQueue::dequeue()
{
    if (isEmpty())
    {
        throw std::underflow_error("Priority queue is empty");
    }
    PriorityNode *temp = front;
    front = front->next;
    delete temp;
}

int PriorityQueue::peek() const
{
    if (isEmpty())
    {
        throw std::underflow_error("Priority queue is empty");
    }
    return front->data;
}

bool PriorityQueue::isEmpty() const
{
    return front == nullptr;
}

void PriorityQueue::display() const
{
    if (isEmpty())
    {
        throw std::underflow_error("Priority queue is empty");
    }
    PriorityNode *current = front;
    std::cout << "Priority Queue: ";
    while (current != nullptr)
    {
        std::cout << current->data << " (Priority: " << current->priority << ")" << " -> ";
        current = current->next;
    }
    std::cout << std::endl;
}