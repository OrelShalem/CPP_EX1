// 313610123
/**
 * @file PriorityQueue.hpp
 * @brief Header file for a Priority Queue implementation
 *
 * This file defines a Priority Queue class that implements a priority queue
 * data structure using a linked list. Elements with higher priority are
 * dequeued before elements with lower priority.
 */

/**
 * @brief Node structure for the priority queue linked list
 *
 * Each node contains the data value, a priority value, and a pointer
 * to the next node in the queue.
 */
struct PriorityNode
{
    int data;           ///< The data stored in the node
    int priority;       ///< The priority of the node (lower value = higher priority)
    PriorityNode *next; ///< Pointer to the next node in the queue
};

/**
 * @brief A class implementing a priority queue data structure
 *
 * This Priority Queue implementation uses a linked list to store elements
 * ordered by their priority. Elements with the highest priority (lowest priority value)
 * are dequeued first.
 */
class PriorityQueue
{
private:
    PriorityNode *front; ///< Pointer to the front of the priority queue

public:
    /**
     * @brief Default constructor initializing an empty priority queue
     */
    PriorityQueue() : front(nullptr) {}

    /**
     * @brief Destructor to free allocated memory
     */
    ~PriorityQueue();

    /**
     * @brief Adds an element with the specified priority to the queue
     *
     * Inserts the element in the appropriate position based on priority.
     *
     * @param value The data value to be added
     * @param priority The priority of the element (lower value = higher priority)
     */
    void enqueue(int value, int priority);

    /**
     * @brief Removes the highest priority element from the queue
     *
     * @throws std::underflow_error If the queue is empty
     */
    void dequeue();

    /**
     * @brief Returns the highest priority element without removing it
     *
     * @return int The highest priority element's value
     * @throws std::underflow_error If the queue is empty
     */
    int peek() const;

    /**
     * @brief Checks if the priority queue is empty
     *
     * @return true If the queue is empty
     * @return false If the queue contains at least one element
     */
    bool isEmpty() const;

    /**
     * @brief Displays all elements in the priority queue
     *
     * Prints all elements along with their priorities to standard output.
     */
    void display() const;
};