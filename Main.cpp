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

void printAdjColorTable(GraphColoring gc, int** matrix)
{
    for (unsigned i = 0; i < gc.nodeNum; i++)
    {
        cout << "line " << i << "\t: ";
        for (unsigned j = 0; j < gc.colorNum; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void printAdjTable(GraphColoring gc, int** matrix)
{
    for (unsigned i = 0; i < gc.nodeNum; i++)
    {
        cout << "line " << i << "\t: ";
        for (unsigned j = 0; j < gc.nodeNum; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void printSolution(GraphColoring& gc,  int* solution)
{
    for (unsigned i = 0; i < gc.nodeNum; i++)
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

int** createAdjList(GraphColoring &gc)
{
    vector<vector<int>> adjList;
    int** adjReturn = new int*[gc.nodeNum];

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

    for (int i = 0; i < gc.nodeNum; i++)
    {
        adjReturn[i] = new int[gc.nodeNum];
    }
    for (int i = 0; i < gc.nodeNum; i++)
    {
        int j = 0;
        for (; j < adjList[i].size(); j++)
        {
            adjReturn[i][j] = adjList[i][j];
        }
        adjReturn[i][j] = -1;
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

    return adjReturn;
}

int* generateInitialSolution(int n, int k)
{
   int* solution = new int[n];
    for (int i = 0; i < n; i++)
    {
        solution[i] = rand(0, k-1);
    }
    return solution;
}


//void solutionToMatrix(GraphColoring &gc, vector<int> solution, vector<vector<int>>& matrix)
//{
//    for (unsigned i = 0; i < gc.edges.size(); i++)
//    {
//        matrix[gc.edges[i][0]][solution[gc.edges[i][1]]]++;
//        matrix[gc.edges[i][1]][solution[gc.edges[i][0]]]++;
//    }
//}

void solutionToMatrix(GraphColoring &gc, int* solution, int** matrix)
{
    // Clean matrix
    for (int i = 0; i < gc.nodeNum; ++i) {
        for (int j = 0; j < gc.colorNum; ++j) {
            matrix[i][j] = 0;
        }
    }


    int a;
    int b;
    for (unsigned i = 0; i < gc.edges.size(); i++)
    {
        a = solution[gc.edges[i][1]];
        b = solution[gc.edges[i][0]];
        matrix[gc.edges[i][0]][solution[gc.edges[i][1]]]++;
        matrix[gc.edges[i][1]][solution[gc.edges[i][0]]]++;
    }
}

int findConflictNode(vector<int> solution, vector<vector<int>> matrix)
{
    int currentConflictValue = 0;
    vector<int> conflictIndex;
    for (unsigned i = 0; i < solution.size(); i++)
    {
        if (matrix[i][solution[i]] > currentConflictValue)
        {
            currentConflictValue = matrix[i][solution[i]];
            conflictIndex.clear();
            conflictIndex.push_back(i);
        }
        else if (matrix[i][solution[i]] == currentConflictValue)
        {
            conflictIndex.push_back(i);
        }
    }
    if (conflictIndex.empty())
    {
        return -1;
    }
    return conflictIndex[rand(0, conflictIndex.size() - 1)];
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

    return choiceIndex[rand(0, choiceIndex.size() - 1)];
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

int calculateDelta(int node, int originalColor, int newColor, int** matrix)
{
    return matrix[node][newColor] - matrix[node][originalColor];
}

int* findMove(GraphColoring &gc, int* solution, int** matrix, int** tabuTable, int iter, int bestDeltaGlobal)
{

    int* returnMove = new int[4];
    int currentMove[3];

    int bestDeltaTabu = gc.nodeNum+1;
    int bestDeltaNonTabu = gc.nodeNum+1;

    int bestDeltaTabuMove[3];
    bestDeltaTabuMove[0] = -1;
    bestDeltaTabuMove[1] = -1;
    bestDeltaTabuMove[2] = -1;

    int bestDeltaNonTabuMove[3];
    bestDeltaNonTabuMove[0] = -1;
    bestDeltaNonTabuMove[1] = -1;
    bestDeltaNonTabuMove[2] = -1;

    int currentDelta;
    int originalColor;
    int currentColor;
    int sampleCount = 1;

    // Every node
    for (int i = 0; i < gc.nodeNum; i++)
    {
        originalColor = solution[i];

        // If there is a conflict Every different color
        // 结构有问题：不需要每一次都遍历
        if (matrix[i][originalColor] > 0) {
            currentColor = -1;

            // 所有颜色
            for (int k = 0; k < gc.colorNum; k++)
            {
                // 和原本颜色一样就跳过
                if (k != originalColor)
                {
                    currentColor = k;
                    currentDelta = matrix[i][currentColor] - matrix[i][originalColor];


                    // Update best if there is
                    // tabu
                    if ((currentDelta <= bestDeltaTabu) && (iter < tabuTable[i][currentColor]))
                    {
                        currentMove[0] = i;
                        currentMove[1] = originalColor;
                        currentMove[2] = currentColor;


                        if (currentDelta < bestDeltaTabu)
                        {
                            bestDeltaTabuMove[0] = currentMove[0];
                            bestDeltaTabuMove[1] = currentMove[1];
                            bestDeltaTabuMove[2] = currentMove[2];

                            bestDeltaTabu = currentDelta;
                            sampleCount = 1;
                        }
                        else if (currentDelta == bestDeltaTabu)
                        {
                            sampleCount++;
                            if (rand(0, sampleCount) == 0)
                            {
                                bestDeltaTabuMove[0] = currentMove[0];
                                bestDeltaTabuMove[1] = currentMove[1];
                                bestDeltaTabuMove[2] = currentMove[2];
                            }
                        }
                    }


                    // non tabu
                    else
                    {
                        currentMove[0] = i;
                        currentMove[1] = originalColor;
                        currentMove[2] = currentColor;

                        if (currentDelta < bestDeltaNonTabu)
                        {
                            bestDeltaNonTabuMove[0] = currentMove[0];
                            bestDeltaNonTabuMove[1] = currentMove[1];
                            bestDeltaNonTabuMove[2] = currentMove[2];
                            bestDeltaNonTabu = currentDelta;
                            sampleCount = 1;
                        }
                        else if (currentDelta == bestDeltaNonTabu)
                        {
                            sampleCount++;
                            if (rand(0, sampleCount) == 0)
                            {
                                bestDeltaNonTabuMove[0] = currentMove[0];
                                bestDeltaNonTabuMove[1] = currentMove[1];
                                bestDeltaNonTabuMove[2] = currentMove[2];
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
        returnMove[0] = bestDeltaTabuMove[0];
        returnMove[1] = bestDeltaTabuMove[1];
        returnMove[2] = bestDeltaTabuMove[2];
        returnMove[3] = bestDeltaTabu;

        return returnMove;
    }

    if (bestDeltaNonTabuMove[0] == -1)
    {
        returnMove[0] = -1;
        returnMove[1] = -1;
        returnMove[2] = -1;
        returnMove[3] = 0;

        return returnMove;
    }

    returnMove[0] = bestDeltaNonTabuMove[0];
    returnMove[1] = bestDeltaNonTabuMove[1];
    returnMove[2] = bestDeltaNonTabuMove[2];
    returnMove[3] = bestDeltaNonTabu;

    return returnMove;
}

void makeMove(GraphColoring& gc, int* move, int &f, int iter, int* solution, int** matrix, int** adj, int** tabuTable, int &bestDeltaGlobal)
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

    for (unsigned i = 0; i < gc.nodeNum; i++)
    {
        if (adj[node][i] == -1)
        {
            break;
        }
        matrix[adj[node][i]][originalColor]--;
        matrix[adj[node][i]][newColor]++;
    }
}

void printMove(int* move)
{
    for (int i = 0; i < 4; i++)
    {
        cout << move[i] << " ";
    }
    cout << endl;
}

int calculateF(GraphColoring &gc, int* solution)
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

void tabuSearch(int maxIter, GraphColoring &gc, int *solution, int** matrix, int** adj)
{
    // Initial iteration = 0
    int iter = 0;
    int f = calculateF(gc, solution);
    int bestDeltaGlobal = gc.nodeNum + 1;
    int* move;

    // Initialize tabu table
    int** tabuTable = new int*[gc.nodeNum];
    for (int i = 0; i < gc.nodeNum; i++)
    {
        tabuTable[i] = new int[gc.colorNum];
    }
    for (int i = 0; i < gc.nodeNum; i++)
    {
        for (int j = 0; j < gc.colorNum; j++)
        {
            tabuTable[i][j] = 0;
        }
    }


    // Until reaches max iteration
    for (iter = 0; iter < maxIter; iter++)
    {
        if (iter % 300000 == 0)
        {
            cout << "iter:\t" << iter << endl;
        }
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
//        cout << iter << "\t";
//        printMove(move);


        if (move[0] == -1)
        {
            return;
        }

        makeMove(gc, move, f, iter, solution, matrix, adj, tabuTable, bestDeltaGlobal);

        if (f == 0)
        {
            cout << iter << endl;
            return;
        }
    }

    // Print information after
//    cout << iter << endl;
//    cout << "matrix: " << endl;
//    printAdjColorTable(gc, matrix);
//    cout << "tabu: " << endl;
//    printAdjColorTable(gc, tabuTable);
//    for (int i = 0; i < gc.nodeNum; i++)
//    {
//        cout << i << ": "
//             << "color: " << solution[i] << "\tconflict: " << matrix[i][solution[i]] << endl;
//    }
//    cout << "f after: " << f << endl;

}

void printColorSet(vector<vector<int>>& colorSet)
{
    for (int i = 0; i < colorSet.size(); i++)
    {
        cout << "Color " << i << "\t:";
        for (int j = 0; j < colorSet[i].size(); j++)
        {
            cout << colorSet[i][j] << " ";
        }
        cout << endl;
    }
}

int findMaxColor(vector<vector<int>>& colorSet)
{
    int setSize = -1;
    int chosenSet = -1;
    int totalSize;
    int currenLineSize;

    totalSize = colorSet.size();
    for (int i = 0; i < totalSize; i++)
    {
        currenLineSize = colorSet[i].size();
        if (currenLineSize > setSize)
        {
            setSize = currenLineSize;
            chosenSet = i;
        }
    }

    return chosenSet;
}

void updateColorSet(int color, int iter, vector<vector<int>>& chosen, vector<vector<int>>& other, vector<vector<int>>& offSpring)
{

    int found = 0;
    int currentNode = -1;

    while (chosen[color].size() > 0)
    {

        // 逐个添加元素
        currentNode = chosen[color].front();
        chosen[color].erase(chosen[color].begin());
        offSpring[iter].push_back(currentNode);

        // 删除另一个的元素
        found = 0;
        for (int j = 0; j < other.size(); j++)
        {
            for (int k = 0; k < other[j].size(); k++)
            {
                if (currentNode == other[j][k])
                {
                    other[j].erase(other[j].begin() + k);
                    found = 1;
                    break;
                }
            }
            if (found == 1)
            {
                break;
            }
        }
    }

}

int* Crossover(GraphColoring& gc, int* solution1, int* solution2)
{
    int maxColor = -1;

    // 创建颜色集
    vector<vector<int>> colorSet1;
    vector<vector<int>> colorSet2;
    vector<vector<int>> offSpringColorSet;

    for (unsigned i = 0; i < gc.colorNum; i++)
    {
        vector<int> line;
        colorSet1.push_back(line);
        colorSet2.push_back(line);
        offSpringColorSet.push_back(line);
    }

    int color1;
    int color2;

    for (int i = 0; i < gc.nodeNum; i++)
    {
        color1 = solution1[i];
        color2 = solution2[i];
        colorSet1[solution1[i]].push_back(i);
        colorSet2[solution2[i]].push_back(i);
    }



    int* offspringSolution = new int[gc.nodeNum];
    for (int i = 0; i < gc.colorNum; i++)
    {
        // 奇数S1
        if (i % 2 == 1)
        {

            maxColor = findMaxColor(colorSet1);
            updateColorSet(maxColor, i, colorSet1, colorSet2, offSpringColorSet);
        }

        // 偶数S2
        else
        {
            maxColor = findMaxColor(colorSet2);
            updateColorSet(maxColor, i, colorSet2, colorSet1, offSpringColorSet);

        }
    }

    // 输入到Solution
    for (int i = 0; i < gc.colorNum; ++i)
    {
        for (int j = 0; j < offSpringColorSet[i].size(); j++)
        {
            offspringSolution[offSpringColorSet[i][j]] = i;
        }
    }

    return offspringSolution;
}

void HEA(int initialPopulationSize, GraphColoring& gc, int** adjList)
{
    int populationSize = initialPopulationSize;
    int bestF = gc.nodeNum * gc.nodeNum;
    int currentF;
    int parent1Index;
    int parent2Index;
    vector<int*> population;
    int* currentSolution;
    int* bestSolution;
    int* offspring;
    int** currentConflictMatrix = new int*[gc.nodeNum];
    for (int i = 0; i < gc.nodeNum; i++)
    {
        currentConflictMatrix[i] = new int[gc.colorNum];
    }


    // Initial population
    for (int i = 0; i < populationSize; ++i)
    {
        population.push_back(generateInitialSolution(gc.nodeNum, gc.colorNum));
    }

    // Initial generation tabu search
    for (int i = 0; i < populationSize; ++i)
    {
        solutionToMatrix(gc, population[i], currentConflictMatrix);
        tabuSearch(5000, gc, population[i], currentConflictMatrix, adjList);

    }

    // Get the best solution from initial population
    for (int i = 0; i < populationSize; ++i) {

        currentF = calculateF(gc, population[i]);

        // 如果已经解决了就直接跳出程序
        if (currentF == 0)
        {
            cout << "fitting solution found." << endl;
            printSolution(gc, population[i]);
            return;
        }

        // 更新最好的S
        if (currentF < bestF)
        {
            bestF = currentF;
        }
    }

    // 开始杂交
    for (int i = 0; i < 100; ++i) {
        cout << "Crossover " << i << endl;

        // Selecting parents, cannot be the same
        parent1Index = rand(0, populationSize);
        parent2Index = rand(0, populationSize);
        while (parent2Index == parent1Index)
        {
            parent2Index = rand(0, populationSize-1);
        }

//        cout << "solution1 " << parent1Index << endl;
//        printSolution(gc, population[parent1Index]);
//        cout << "solution2" << parent2Index << endl;
//        printSolution(gc, population[parent2Index]);
        cout << populationSize << endl;
        cout << parent1Index << "\t" << parent2Index << endl;

        //crossover
        offspring = Crossover(gc, population[parent1Index], population[parent2Index]);
        solutionToMatrix(gc, offspring, currentConflictMatrix);
        tabuSearch(5000, gc, offspring, currentConflictMatrix, adjList);

        // Update populations
        currentF = calculateF(gc, offspring);
        if (currentF == 0)
        {
            cout << "fitting solution found." << endl;
            printSolution(gc, offspring);
            return;
        }

        // 更新最好的S
        if (currentF < bestF)
        {
            bestF = currentF;
            bestSolution = population[i];
        }

        population.push_back(offspring);
        populationSize++;

    }
    cout << bestF << endl;


}










int main(int argc, char *argv[])
{

    clock_t start, end;
    start = clock();


    // Read input directly from argv[1] and argv[2]
    long long secTimeout = atoll(argv[1]);
    int randSeed = atoi(argv[2]);
    // cout << secTimeout << endl;
    // cout << randSeed << endl;


    // Random Number
    initRand(randSeed);


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

    int** adjList = createAdjList(gc);


    // Creating 2d vector for Adjacent Color Table
    int** adjacentColorTable = new int*[gc.nodeNum];
    for (int i = 0; i < gc.nodeNum; i++)
    {
        adjacentColorTable[i] = new int[gc.colorNum];
    }
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
    int *solution = new int[gc.nodeNum];
    for (int i = 0; i < gc.nodeNum; i++)
    {
        solution[i] = rand(0, gc.colorNum);
    }

    int *solution2 = new int[gc.nodeNum];
    for (int i = 0; i < gc.nodeNum; i++)
    {
        solution2[i] = rand(0, gc.colorNum);
    }

    /*
     * P是什么
     * 禁忌的迭代次数应该是多少次
     * 池应该怎么更新，好像就是直接加进去，那没事了
     */

//     cout << gc.nodeNum << endl;
//     for (unsigned i = 0; i < gc.nodeNum; i++)
//     {
//         cout << solution[i] << endl;
//     }


    // Solution to matrix
    solutionToMatrix(gc, solution, adjacentColorTable);
//     for (unsigned i = 0; i < gc.nodeNum; i++)
//     {
//         cout << "line " << i << ": ";
//         for (unsigned j = 0; j < gc.colorNum; j++)
//         {
//             cout << adjacentColorTable[i][j] << " ";
//         }
//         cout << endl;
//     }



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

//    tabuSearch(9999999, gc, solution, adjacentColorTable, adjList);
    HEA(5, gc, adjList);

    end = clock();
    cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

    return 0;
}