#include "UnionFind.hpp"
#include <stdexcept>

UnionFind::UnionFind(int n) : parent(new int[n]), size(new int[n]), setCount(n)
{
    for (int i = 0; i < n; i++)
    {
        parent[i] = i;
        size[i] = 1;
    }
}

UnionFind::~UnionFind()
{
    delete[] parent;
    delete[] size;
}

int UnionFind::find(int x)
{
    if (parent[x] != x)
    {
        parent[x] = find(parent[x]);
    }
    return parent[x];
}

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

bool UnionFind::connected(int x, int y)
{
    return find(x) == find(y);
}

int UnionFind::getSetCount() const
{
    return setCount;
}
