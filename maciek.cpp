#include<iostream>
#include<vector>
#include<map>
#include<set>

using namespace std;

struct Edge {
    int v1;
    int v2;
    double weight;
};

struct Node {
    int index;
    vector<int> edgesIndices;
    double label;
    bool visited;
};

void initializeLabels(vector<Node>& nodes, vector<Edge>& edges) {
    dfsInitializeLabels(nodes, edges, 0, true);
}

void dfsInitializeLabels(vector<Node>& nodes, vector<Edge>& edges, int nodeIndex, bool isHouse) {
    if (nodes[nodeIndex].visited)
        return;
    nodes[nodeIndex].visited = true;

    if (isHouse)
        nodes[nodeIndex].label = 0.0;
    else {
        double maxNeighborWeight = 0.0;
        for (auto i: nodes[nodeIndex].edgesIndices)
            maxNeighborWeight = max(maxNeighborWeight, edges[i].weight);
        nodes[nodeIndex].label = maxNeighborWeight;
    }

    for (auto i: nodes[nodeIndex].edgesIndices) {
        dfsInitializeLabels(nodes, edges, edges[i].v1 == nodeIndex ? edges[i].v2 : edges[i].v1, !isHouse);
    }
}

vector<Node> constructEqualityGraph(const vector<Node> nodes, const vector<Edge> edges) {
    map<int, Node> equalityGraphMap;

    for (int nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++) {
        for (int edgeIndex: nodes[nodeIndex].edgesIndices) {
            Edge edge = edges[edgeIndex];
            if (nodes[edge.v1].label + nodes[edge.v2].label == edge.weight) {
                if (equalityGraphMap.find(nodeIndex) != equalityGraphMap.end()) {
                    Node newNode;
                    newNode.index = nodeIndex;
                    newNode.edgesIndices = { edgeIndex };
                    newNode.label = nodes[nodeIndex].label;
                    newNode.visited = false;
                    equalityGraphMap[nodeIndex] = newNode;
                } else {
                    equalityGraphMap[edge.v1].edgesIndices.push_back(edgeIndex);
                }
            }
        }
    }

    vector<Node> result;
    for (const auto& pair : equalityGraphMap) {
        result.push_back(pair.second);
    }

    return result;
}

void resetVisited(vector<Node>& nodes) {
    for (auto i : nodes) {
        i.visited = false;
    }
}

int findNotAssignedWell(vector<Node>& nodes, const vector<Edge>& edges, const vector<Node>& equalityGraph) {
    resetVisited(nodes);
    set<int> verticesInEqualityGraph;
    for (auto i: equalityGraph) {
        verticesInEqualityGraph.insert(i.index);
    }

    return dfsFindNotAssignedWell(nodes, edges, verticesInEqualityGraph, 0, false);
}

int dfsFindNotAssignedWell(vector<Node>& nodes, const vector<Edge>& edges, const set<int> verticesInEqualityGraph, int nodeIndex, int isWell) {
    if (nodes[nodeIndex].visited)
        return -1;
    nodes[nodeIndex].visited = true;
    // Tu pewnie nie musi zwracać studni ale tak łatwiej się o tym myśli
    if (isWell && verticesInEqualityGraph.find(nodeIndex) != verticesInEqualityGraph.end())
        return nodeIndex;
    
    int result = -1;
    for (auto i : nodes[nodeIndex].edgesIndices) {
        result = dfsFindNotAssignedWell(nodes, edges, verticesInEqualityGraph, edges[i].v1 == nodeIndex ? edges[i].v2 : edges[i].v1, !isWell);
        if (result != -1)
            break;
    }
    return result;
}

int findIndexOfVertex(const vector<Node>& nodes, int vertexIndex) {
    for (int i = 0; i < nodes.size(); i++) {
        if (nodes[i].index == vertexIndex)
            return i;
    }
    return -1;
}

vector<vector<int>> hunigerianAlgorithm(vector<Node> nodes, vector<Edge> edges) {
    resetVisited(nodes);
    initializeLabels(nodes, edges);
    vector<Node> equalityGraph = constructEqualityGraph(nodes, edges);
    
    vector<Node> currentMatching;
    while(currentMatching.size() < nodes.size()/2) {
        int alternatingPathStart = findNotAssignedWell(nodes, edges, equalityGraph);
        vector<Node> alternatingPath;
        int indexInEqualityGraph = findIndexOfVertex(equalityGraph, alternatingPathStart);
        Node startNode = equalityGraph[indexInEqualityGraph];
        alternatingPath.push_back(startNode);
        while(alternatingPath.size() % 2 == 0) {
            
        }
    }

    return {};
}

int main() {
    cout << "Hello world" << endl;
}