/**
 * @file UnionFind.hpp
 * @brief Header file for a Union-Find (Disjoint-Set) data structure
 *
 * This file defines a Union-Find class which implements the disjoint-set
 * data structure. It supports operations to find which set an element belongs to,
 * unite two sets, and check if two elements are in the same set.
 */

/**
 * @brief A class implementing the Union-Find (Disjoint-Set) data structure
 *
 * This implementation uses path compression in the find operation and
 * union by size for the merge operation to achieve near-constant time operations.
 */
class UnionFind
{
private:
    int *parent;  ///< Array storing parent pointers for each element
    int *size;    ///< Array storing the size of each set
    int setCount; ///< Number of disjoint sets

public:
    /**
     * @brief Constructs a Union-Find data structure with n elements
     *
     * Initializes n singleton sets, each containing one element.
     *
     * @param n Number of elements
     */
    UnionFind(int n);

    /**
     * @brief Destructor to free allocated memory
     */
    ~UnionFind();

    /**
     * @brief Finds the representative (root) of the set containing x
     *
     * Uses path compression to optimize future queries.
     *
     * @param x The element to find
     * @return int The representative of the set containing x
     */
    int find(int x);

    /**
     * @brief Unites the sets containing elements x and y
     *
     * Uses union by size to keep the tree balanced.
     *
     * @param x First element
     * @param y Second element
     */
    void unionSet(int x, int y);

    /**
     * @brief Checks if elements x and y are in the same set
     *
     * @param x First element
     * @param y Second element
     * @return true If x and y are in the same set
     * @return false If x and y are in different sets
     */
    bool connected(int x, int y);

    /**
     * @brief Gets the current number of disjoint sets
     *
     * @return int Number of disjoint sets
     */
    int getSetCount() const;
};