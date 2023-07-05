#include <iostream>
#include <string>
#include <chrono>
#include <array>
#include <queue>
#include <vector>
#include <functional>
#include <random>

using namespace std;
using NodeId = int;
using EdgeId = NodeId;
using ColorId = NodeId;

using Edge = std::array<NodeId, 2>; // undirected link.
using move = std::array<int, 3>;
using line = array<int, 5>;

struct GraphColoring
{
    NodeId nodeNum;
    EdgeId edgeNum;
    ColorId colorNum;
    std::vector<Edge> edges;
};

mt19937 pseudoRandNumGen;
void initRand(int seed) { pseudoRandNumGen = mt19937(seed); }
int fastRand(int lb, int ub) { return (pseudoRandNumGen() % (ub - lb)) + lb; }
int fastRand(int ub) { return pseudoRandNumGen() % ub; }
int rand(int lb, int ub) { return uniform_int_distribution<int>(lb, ub - 1)(pseudoRandNumGen); }
int rand(int ub) { return uniform_int_distribution<int>(0, ub - 1)(pseudoRandNumGen); }

using NodeColors = std::vector<ColorId>;

int generateRandomInt(int lowerBound, int upperBound)
{
    srand(time(nullptr) + clock());
    int randNum = (rand() % (upperBound - lowerBound + 1)) + lowerBound;
    return randNum;
}

void printAdjColorTable(vector<vector<int>> matrix)
{
    for (unsigned i = 0; i < matrix.size(); i++)
    {
        cout << "line " << i << "\t: ";
        for (unsigned j = 0; j < matrix[i].size(); j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void printSolution(vector<int> solution)
{
    for (unsigned i = 0; i < solution.size(); i++)
    {
        cout << "line " << i << "\t: " << solution[i] << endl;
    }
}

void loadInput(istream &is, GraphColoring &gc)
{
    is >> gc.nodeNum >> gc.edgeNum >> gc.colorNum;
    gc.edges.resize(gc.edgeNum);
    for (auto edge = gc.edges.begin(); edge != gc.edges.end(); ++edge)
    {
        is >> (*edge)[0] >> (*edge)[1];
    }
}

void saveOutput(ostream &os, NodeColors &nodeColors)
{
    for (auto color = nodeColors.begin(); color != nodeColors.end(); ++color)
    {
        os << *color << endl;
    }
}

vector<vector<int>> createACT(int n, int k)
{
    vector<vector<int>> matrix;
    for (unsigned i = 0; i < n; i++)
    {
        vector<int> innerMatrix(k, 0);
        matrix.push_back(innerMatrix);
    }
    return matrix;
}

vector<vector<int>> createAdjList(GraphColoring &gc)
{
    vector<vector<int>> adjList;
    for (unsigned i = 0; i < gc.nodeNum; i++)
    {
        vector<int> line;
        adjList.push_back(line);
    }
    for (unsigned i = 0; i < gc.edges.size(); i++)
    {
        // cout << gc.edges[i][0] << " " << gc.edges[i][1] << endl;
        adjList[gc.edges[i][0]].push_back(gc.edges[i][1]);
        adjList[gc.edges[i][1]].push_back(gc.edges[i][0]);
    }

    // cout << adjList.size() << endl;
    // for (int i = 0; i < adjList.size(); i++)
    // {
    //     for (int j = 0; j < adjList[i].size(); j++)
    //     {
    //         cout << adjList[i][j] << "\t";
    //     }
    //     cout << endl;
    // }

    return adjList;
}

vector<int> generateInitialSolution(int n, int k)
{
    vector<int> solution(n, 0);
    for (int i = 0; i < solution.size(); i++)
    {
        solution[i] = generateRandomInt(0, k - 1);
    }
    return solution;
}

void solutionToMatrix(GraphColoring &gc, vector<int> solution, vector<vector<int>> &matrix)
{
    for (unsigned i = 0; i < gc.edges.size(); i++)
    {
        matrix[gc.edges[i][0]][solution[gc.edges[i][1]]]++;
        matrix[gc.edges[i][1]][solution[gc.edges[i][0]]]++;
    }
}

int findConflictNode(vector<int> soluton, vector<vector<int>> matrix)
{
    int currentConflictValue = 0;
    vector<int> conflictIndex;
    for (unsigned i = 0; i < soluton.size(); i++)
    {
        if (matrix[i][soluton[i]] > currentConflictValue)
        {
            currentConflictValue = matrix[i][soluton[i]];
            conflictIndex.clear();
            conflictIndex.push_back(i);
        }
        else if (matrix[i][soluton[i]] == currentConflictValue)
        {
            conflictIndex.push_back(i);
        }
    }
    if (conflictIndex.empty())
    {
        return -1;
    }
    return conflictIndex[generateRandomInt(0, conflictIndex.size() - 1)];
}

int findNewColor(vector<int> solution, vector<vector<int>> matrix, int nodeIndex)
{
    int oldColor = solution[nodeIndex];
    vector<int> choiceIndex;
    int currentChoice = 999999;
    for (unsigned i = 1; i < matrix[nodeIndex].size(); i++)
    {
        if (i == oldColor)
        {
            continue;
        }

        if (matrix[nodeIndex][i] < currentChoice)
        {
            currentChoice = matrix[nodeIndex][i];
            choiceIndex.clear();
            choiceIndex.push_back(i);
        }
        else if (matrix[nodeIndex][i] == currentChoice)
        {
            choiceIndex.push_back(i);
        }
    }

    return choiceIndex[generateRandomInt(0, choiceIndex.size() - 1)];
}

void updateSolAndMatrix(vector<int> &solution, vector<vector<int>> &matrix, vector<vector<int>> adj, int node, int choice)
{
    int originalColor = solution[node];
    int newColor = choice;

    solution[node] = choice;

    for (unsigned i = 0; i < adj[node].size(); i++)
    {
        if (adj[node][i] == 1)
        {
            matrix[i][originalColor]--;
            matrix[i][newColor]++;
        }
    }
}

vector<vector<int>> createTabuTable(GraphColoring &gc)
{
    vector<vector<int>> tabuTable;
    for (unsigned i = 0; i < gc.nodeNum; i++)
    {
        vector<int> line(gc.colorNum, 0);
        tabuTable.push_back(line);
    }
    return tabuTable;
}

int calculateDelt(int node, int originalColor, int newColor, vector<vector<int>> matrix)
{
    return matrix[node][newColor] - matrix[node][originalColor];
}

array<int, 4> findMove(GraphColoring &gc, vector<int> solution, vector<vector<int>> matrix, vector<vector<int>> tabuTable, int iter, int bestDeltaGlobal)
{
    using move = std::array<int, 3>;
    std::array<int, 3> currentMove{};
    std::array<int, 4> returnMove{};

    int bestDeltaTabu = gc.nodeNum+1;
    int bestDeltaNonTabu = gc.nodeNum+1;
    std::array<int, 3> bestDeltaTabuMove = {-1, -1, -1};
    std::array<int, 3> bestDeltaNonTabuMove = {-1, -1, -1};
    int currentDelta;
    int originalColor;
    int currentColor;

    // Every node
    for (int i = 0; i < solution.size(); i++)
    {
        originalColor = solution[i];

        // If there is a conflict Every different color
        // 结构有问题：不需要每一次都遍历
        if (matrix[i][originalColor] > 0) {
            currentColor = -1;

            // 所有颜色
            for (int k = 0; k < gc.colorNum; k++) {
                // 和原本颜色一样就跳过
                if (k != originalColor)
                {

                    currentDelta = calculateDelt(i, originalColor, k, matrix);
                    currentColor = k;
                    currentMove = {i, originalColor, currentColor};

                    // Update best if there is
                    // tabu
                    if (iter < tabuTable[i][currentColor]) {
                        if (currentDelta < bestDeltaTabu) {
                            bestDeltaTabuMove = currentMove;
                            bestDeltaTabu = currentDelta;
                        } else if (currentDelta == bestDeltaTabu) {
                            if (generateRandomInt(0, 1) == 1) {
                                bestDeltaTabu = currentDelta;
                            }
                        }
                    }

                    // non tabu
                    else
                    {
                        if (currentDelta < bestDeltaNonTabu)
                        {
                            bestDeltaNonTabuMove = currentMove;
                            bestDeltaNonTabu = currentDelta;
                        } else if (currentDelta == bestDeltaNonTabu)
                        {
                            if (generateRandomInt(0, 1) == 1)
                            {
                                bestDeltaNonTabu = currentDelta;
                            }
                        }

                    }
                }
            }
        }
    }

    // Which to return
    if (bestDeltaTabu < bestDeltaGlobal && bestDeltaTabu < bestDeltaNonTabu)
    {
        if (bestDeltaTabuMove[0] == -1)
        {
            returnMove = {-1, -1, -1, 0};
            return returnMove;
        }

        returnMove = {bestDeltaTabuMove[0], bestDeltaTabuMove[1], bestDeltaTabuMove[2], bestDeltaTabu};
        return returnMove;
    }

    if (bestDeltaNonTabuMove[0] == -1)
    {
        returnMove = {-1, -1, -1, 0};
        return returnMove;
    }

    returnMove = {bestDeltaNonTabuMove[0], bestDeltaNonTabuMove[1], bestDeltaNonTabuMove[2], bestDeltaNonTabu};
    return returnMove;
}

void makeMove(array<int, 4> move, int &f, int iter, vector<int> &solution, vector<vector<int>> &matrix, vector<vector<int>> adj, vector<vector<int>> &tabuTable, int &bestDeltaGlobal)
{
    int node = move[0];
    int originalColor = move[1];
    int newColor = move[2];
    int delta = move[3];

    if (delta < bestDeltaGlobal)
    {
        bestDeltaGlobal = delta;
    }
    solution[node] = newColor;
    f = f + delta;
    tabuTable[node][originalColor] = iter + f + rand() % 10;

    for (unsigned i = 0; i < adj[node].size(); i++)
    {

        matrix[adj[node][i]][originalColor]--;
        matrix[adj[node][i]][newColor]++;
    }
}


void printMove(array<int, 4> move)
{
    for (int i = 0; i < move.size(); i++)
    {
        cout << move[i] << " ";
    }
    cout << endl;
}

int calculateInitialF(GraphColoring &gc, vector<int> solution)
{
    int f = 0;
    for (unsigned i = 0; i < gc.edgeNum; i++)
    {
        if (solution[gc.edges[i][0]] == solution[gc.edges[i][1]])
        {
            f++;
        }
    }

    return f;
}

void tabuSearch(GraphColoring &gc, vector<int> &solution, vector<vector<int>> &matrix, vector<vector<int>> adj)
{
    // Initial iteration = 0
    int iter = 0;
    int f = calculateInitialF(gc, solution);
    int bestDeltaGlobal = gc.nodeNum + 1;
    array<int, 4> move;

    // Initialize tabu table
    vector<vector<int>> tabuTable = createTabuTable(gc);

    for (iter = 0; iter < 10000; iter++)
    {
        move = findMove(gc, solution, matrix, tabuTable, iter, bestDeltaGlobal);

        // cout << "\n\n\n\n\niter = " << iter << endl;
        // cout << "matrix: " << endl;
        // printAdjColorTable(matrix);
        // cout << "tabu: " << endl;
        // printAdjColorTable(tabuTable);
        // for (int i = 0; i < solution.size(); i++)
        // {
        //     cout << i << ": " << "color: " << solution[i] << "\tconflct: " << matrix[i][solution[i]] << endl;
        // }

        // cout << "f before: " << f << endl;
        // cout << "bestDeltaGlobal: " << bestDeltaGlobal << endl;
        cout << iter << "\t";
        printMove(move);


        if (move[0] == -1)
        {
            break;
        }

        makeMove(move, f, iter, solution, matrix, adj, tabuTable, bestDeltaGlobal);

        if (f == 0)
        {
            break;
        }
    }

    // Print information after
    cout << iter << endl;
    cout << "matrix: " << endl;
    printAdjColorTable(matrix);
    cout << "tabu: " << endl;
    printAdjColorTable(tabuTable);
    for (int i = 0; i < solution.size(); i++)
    {
        cout << i << ": "
             << "color: " << solution[i] << "\tconflict: " << matrix[i][solution[i]] << endl;
    }
    cout << "f after: " << f << endl;
}

int main(int argc, char *argv[])
{

    clock_t start, end;
    start = clock();

    mt19937(seed);

    // Read input directly from argv[1] and argv[2]
    long long secTimeout = atoll(argv[1]);
    int randSeed = atoi(argv[2]);
    // cout << secTimeout << endl;
    // cout << randSeed << endl;

    // Creating struct and read input
    GraphColoring gc;
    loadInput(cin, gc);
    // cout << gc.nodeNum << endl;
    // cout << gc.edgeNum << endl;
    // cout << gc.colorNum << endl;
    // cout << gc.edges.size() << endl;
    // for (unsigned i = 0; i < gc.edges.size(); i++)
    // {
    //     cout << gc.edges[i][0] << " " << gc.edges[i][1] << endl;
    // }

    // Generating adjList for future search
    vector<vector<int>> adjList = createAdjList(gc);

    // Creating 2d vector for Adjacent Color Table
    vector<vector<int>> adjacentColorTable = createACT(gc.nodeNum, gc.colorNum);
    // cout << adjacentColorTable.size() << endl;
    // for (unsigned i = 0; i < adjacentColorTable.size(); i++)
    // {
    //     cout << "line " << i << ": ";
    //     for (unsigned j = 0; j < adjacentColorTable[i].size(); j++)
    //     {
    //         cout << adjacentColorTable[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // Generating initial solution
    vector<int> solution = generateInitialSolution(gc.nodeNum, gc.colorNum);
    // cout << solution.size() << endl;
    // for (unsigned i = 0; i < solution.size(); i++)
    // {
    //     cout << solution[i] << endl;
    // }

    // Solution to matrix
    solutionToMatrix(gc, solution, adjacentColorTable);
    // for (unsigned i = 0; i < adjacentColorTable.size(); i++)
    // {
    //     cout << "line " << i << ": ";
    //     for (unsigned j = 0; j < adjacentColorTable[i].size(); j++)
    //     {
    //         cout << adjacentColorTable[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // Find current conflict
    // int conflictIndex;
    // int newColor;
    // printAdjColorTable(adjacentColorTable);
    // printSolution(solution);

    // for (size_t i = 0; i < 10000; i++)
    // {
    //     conflictIndex = findConflictNode(solution, adjacentColorTable);
    //     if (conflictIndex == -1)
    //     {
    //         break;
    //     }
    //     newColor = findNewColor(solution, adjacentColorTable, conflictIndex);
    //     cout << "Node: " << conflictIndex << " Original Color: " << solution[conflictIndex] << " New Color: " << newColor << endl;
    //     for (unsigned i = 0; i < adjacentColorTable[conflictIndex].size(); i++)
    //     {
    //         cout << adjacentColorTable[conflictIndex][i] << " ";
    //     }
    //     cout << endl;
    //     updateSolAndMatrix(solution, adjacentColorTable, adjList, conflictIndex, newColor);
    //     printAdjColorTable(adjacentColorTable);
    // }

    // cout << "Conflicting color from neighbour: " << endl;
    // for (size_t i = 0; i < solution.size(); i++)
    // {
    //     cout << "line " << i << ": ";
    //     cout << adjacentColorTable[i][solution[i]] << endl;
    // }

    tabuSearch(gc, solution, adjacentColorTable, adjList);

    end = clock();
    cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

    return 0;
}