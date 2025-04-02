/**
 * @file UnionFind.cpp
 * @brief Implementation of the UnionFind class
 */
#include "UnionFind.hpp"
#include <stdexcept>

/**
 * @brief Constructs a Union-Find data structure with n elements
 *
 * Initializes n singleton sets, where each element is its own parent
 * and each set has a size of 1.
 *
 * @param n Number of elements
 */
UnionFind::UnionFind(int n) : parent(new int[n]), size(new int[n]), setCount(n)
{
    for (int i = 0; i < n; i++)
    {
        parent[i] = i;
        size[i] = 1;
    }
}

/**
 * @brief Destructor to free allocated memory
 *
 * Deallocates the parent and size arrays.
 */
UnionFind::~UnionFind()
{
    delete[] parent;
    delete[] size;
}

/**
 * @brief Finds the representative (root) of the set containing x
 *
 * Uses path compression to optimize future queries by making every
 * visited node point directly to the root.
 *
 * @param x The element to find
 * @return int The representative of the set containing x
 */
int UnionFind::find(int x)
{
    if (parent[x] != x)
    {
        parent[x] = find(parent[x]);
    }
    return parent[x];
}

/**
 * @brief Unites the sets containing elements x and y
 *
 * If x and y are already in the same set, does nothing.
 * Otherwise, merges the smaller set into the larger one
 * (union by size) and decrements the set count.
 *
 * @param x First element
 * @param y Second element
 */
void UnionFind::unionSet(int x, int y)
{
    int rootX = find(x);
    int rootY = find(y);

    if (rootX == rootY)
    {
        return;
    }

    if (size[rootX] < size[rootY])
    {
        parent[rootX] = rootY;
        size[rootY] += size[rootX];
    }
    else
    {
        parent[rootY] = rootX;
        size[rootX] += size[rootY];
    }
    setCount--;
}

/**
 * @brief Checks if elements x and y are in the same set
 *
 * Two elements are in the same set if they have the same representative.
 *
 * @param x First element
 * @param y Second element
 * @return true If x and y are in the same set
 * @return false If x and y are in different sets
 */
bool UnionFind::connected(int x, int y)
{
    return find(x) == find(y);
}

/**
 * @brief Gets the current number of disjoint sets
 *
 * @return int Number of disjoint sets
 */
int UnionFind::getSetCount() const
{
    return setCount;
}
