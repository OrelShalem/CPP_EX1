// 313610123
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
    // Create a new node and initialize its fields
    PriorityNode *newNode = new PriorityNode;
    newNode->data = value;
    newNode->priority = priority;
    newNode->next = nullptr;

    // If queue is empty or new node has higher priority than front node
    if (isEmpty() || priority < front->priority)
    {
        // Insert at the beginning of the queue
        newNode->next = front;
        front = newNode;
    }
    else
    {
        // Find the position to insert the new node
        PriorityNode *current = front;
        // Traverse until we find a node with lower priority (higher value)
        // or reach the end of the queue
        while (current->next != nullptr && current->next->priority <= priority)
        {
            current = current->next;
        }
        // Insert the new node after the current node
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
    // Check if the queue is empty
    if (isEmpty())
    {
        throw std::underflow_error("Priority queue is empty");
    }
    // Save the front node to be deleted
    PriorityNode *temp = front;
    // Update the front pointer to the next node
    front = front->next;
    // Free the memory of the removed node
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
    // Check if the queue is empty
    if (isEmpty())
    {
        throw std::underflow_error("Priority queue is empty");
    }
    // Return the data of the front node (highest priority)
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
    // Queue is empty if front pointer is null
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
    // Check if the queue is empty
    if (isEmpty())
    {
        throw std::underflow_error("Priority queue is empty");
    }
    // Start from the front node
    PriorityNode *current = front;
    std::cout << "Priority Queue: ";
    // Traverse the queue and print each node's data and priority
    while (current != nullptr)
    {
        std::cout << current->data << " (Priority: " << current->priority << ")" << " -> ";
        current = current->next;
    }
    std::cout << std::endl;
}