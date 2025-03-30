class UnionFind
{
private:
    int *parent;
    int *size;
    int setCount;

public:
    UnionFind(int n);
    ~UnionFind();
    int find(int x);
    void unionSet(int x, int y);
    bool connected(int x, int y);
    int getSetCount() const;
};