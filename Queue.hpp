// orel8155@gmail.com
/**
 * @file Queue.hpp
 * @brief Header file for a simple Queue implementation
 *
 * This file defines a Queue class that implements a basic queue data structure
 * using an array. It provides the standard queue operations like enqueue, dequeue,
 * peek, and status checks.
 */

#include <iostream>
using namespace std;

/**
 * @brief A class implementing a basic queue data structure
 *
 * This Queue implementation uses a circular array to store elements
 * and provides standard queue operations with constant time complexity.
 */
class Queue
{
private:
    int *arr;     ///< Array to store queue elements
    int front;    ///< Index of the front element
    int rear;     ///< Index of the rear element
    int capacity; ///< Maximum capacity of the queue
    int count;    ///< Current number of elements in the queue

public:
    /**
     * @brief Constructs a queue with the specified capacity
     *
     * @param size Maximum number of elements the queue can hold
     */
    Queue(int size);

    /**
     * @brief Destructor to free allocated memory
     */
    ~Queue();

    /**
     * @brief Adds an element to the back of the queue
     *
     * @param x Element to be added
     * @throws std::overflow_error If the queue is full
     */
    void enqueue(int x);

    /**
     * @brief Removes and returns the element at the front of the queue
     *
     * @return int The front element
     * @throws std::underflow_error If the queue is empty
     */
    int dequeue();

    /**
     * @brief Returns the element at the front without removing it
     *
     * @return int The front element
     * @throws std::underflow_error If the queue is empty
     */
    int peek();

    /**
     * @brief Gets the current number of elements in the queue
     *
     * @return int Current queue size
     */
    int getSize();

    /**
     * @brief Checks if the queue is empty
     *
     * @return true If the queue is empty
     * @return false If the queue contains at least one element
     */
    bool isEmpty();

    /**
     * @brief Checks if the queue is full
     *
     * @return true If the queue is full
     * @return false If the queue has space for at least one more element
     */
    bool isFull();
};
