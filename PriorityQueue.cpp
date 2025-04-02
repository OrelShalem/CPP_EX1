/**
 * @file PriorityQueue.cpp
 * @brief Implementation of the PriorityQueue class
 */
#include "PriorityQueue.hpp"
#include <stdexcept>
#include <iostream>

/**
 * @brief Destructor to free allocated memory
 *
 * Deallocates all nodes in the priority queue by repeatedly
 * calling dequeue until the queue is empty.
 */
PriorityQueue::~PriorityQueue()
{
    while (!isEmpty())
    {
        dequeue();
    }
}

/**
 * @brief Adds an element with the specified priority to the queue
 *
 * Creates a new node with the given value and priority, and inserts it
 * into the queue in the correct position based on priority (lower value = higher priority).
 * Elements with the same priority are arranged in FIFO order.
 *
 * @param value The data value to be added
 * @param priority The priority of the element
 */
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

/**
 * @brief Removes the highest priority element from the queue
 *
 * Removes the element at the front of the queue (highest priority)
 * and frees the memory allocated for its node.
 *
 * @throws std::underflow_error If the queue is empty
 */
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

/**
 * @brief Returns the highest priority element without removing it
 *
 * @return int The value of the highest priority element
 * @throws std::underflow_error If the queue is empty
 */
int PriorityQueue::peek() const
{
    if (isEmpty())
    {
        throw std::underflow_error("Priority queue is empty");
    }
    return front->data;
}

/**
 * @brief Checks if the priority queue is empty
 *
 * @return true If the queue is empty
 * @return false If the queue contains at least one element
 */
bool PriorityQueue::isEmpty() const
{
    return front == nullptr;
}

/**
 * @brief Displays all elements in the priority queue
 *
 * Prints each element along with its priority to the standard output
 * in the format: "value (Priority: priority) -> ".
 *
 * @throws std::underflow_error If the queue is empty
 */
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