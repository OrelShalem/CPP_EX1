// orel8155@gmail.com
/**
 * @file Queue.cpp
 * @brief Implementation of the Queue class
 */
#include "Queue.hpp"
#include <iostream>

/**
 * @brief Constructs a queue with the specified capacity
 *
 * Initializes a new queue with the given size, allocating
 * memory for the internal array.
 *
 * @param size Maximum number of elements the queue can hold
 */
Queue::Queue(int size) : arr(new int[size]), front(0), rear(-1), capacity(size), count(0) {}

/**
 * @brief Destructor to free allocated memory
 *
 * Frees the memory allocated for the internal array.
 */
Queue::~Queue()
{
    delete[] arr;
}

/**
 * @brief Adds an element to the back of the queue
 *
 * Adds the specified item to the rear of the queue if space is available.
 * Uses a circular array implementation to efficiently wrap around.
 *
 * @param item Element to be added
 * @throws std::overflow_error If the queue is full
 */
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

/**
 * @brief Removes and returns the element at the front of the queue
 *
 * Removes the element at the front of the queue and returns it.
 * Advances the front index using a circular array implementation.
 *
 * @return int The front element
 * @throws std::underflow_error If the queue is empty
 */
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

/**
 * @brief Returns the element at the front without removing it
 *
 * Provides access to the front element without modifying the queue.
 *
 * @return int The front element
 * @throws std::underflow_error If the queue is empty
 */
int Queue::peek()
{
    if (isEmpty())
    {
        throw std::underflow_error("Queue is empty");
    }
    return arr[front];
}

/**
 * @brief Checks if the queue is empty
 *
 * @return true If the queue is empty
 * @return false If the queue contains at least one element
 */
bool Queue::isEmpty()
{
    return count == 0;
}

/**
 * @brief Checks if the queue is full
 *
 * @return true If the queue is full
 * @return false If the queue has space for at least one more element
 */
bool Queue::isFull()
{
    return count == capacity;
}

/**
 * @brief Gets the current number of elements in the queue
 *
 * @return int Current queue size
 */
int Queue::getSize()
{
    return count;
}