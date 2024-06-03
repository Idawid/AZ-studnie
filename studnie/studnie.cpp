#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <cassert>
#include <algorithm>
#include <filesystem>

#define E 10000

const int INF = 1e9;

using namespace std;
namespace fs = std::filesystem;

struct Point {
    double x, y;
};

double randomDouble(double min, double max) {
    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<double> dis(min, max);
    return int(dis(gen));
}

double euclideanDistance(Point p1, Point p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// Function to generate random input and write it to a file
void generateRandomInput(const string& filename) {
    // Random number of wells (n) and houses each well supplies (k)
    int n = randomDouble(1, 5); // Change range as needed
    int k = randomDouble(1, 5);  // Change range as needed

    // Adjust n and k so that n * k is <= 10
    while (n * k > 10) {
        n = randomDouble(1, 5); // Change range as needed
        k = randomDouble(1, 5);  // Change range as needed
    }

    // Random coordinates for wells
    vector<Point> wells(n);
    for (int i = 0; i < n; ++i) {
        wells[i] = { randomDouble(0, 10), randomDouble(0, 10) }; // Change range as needed
    }

    // Random coordinates for houses
    vector<Point> houses(n * k);
    for (int i = 0; i < n * k; ++i) {
        houses[i] = { randomDouble(0, 10), randomDouble(0, 10) }; // Change range as needed
    }

    // Write random input to the file
    ofstream file(filename);
    if (file.is_open()) {
        file << n << " " << k << endl;
        for (const auto& well : wells) {
            file << well.x << " " << well.y << endl;
        }
        for (const auto& house : houses) {
            file << house.x << " " << house.y << endl;
        }
        file.close();
        cout << "Random input generated and written to " << filename << endl;
    }
    else {
        cerr << "Error writing random input to file: " << filename << endl;
    }
}

// Function to find the minimum cost in the cost matrix
double findMinCost(const vector<vector<double>>& costMatrix, vector<bool>& rowCovered, vector<bool>& colCovered) {
    double minCost = numeric_limits<double>::infinity();
    for (int i = 0; i < costMatrix.size(); ++i) {
        if (!rowCovered[i]) {
            for (int j = 0; j < costMatrix[i].size(); ++j) {
                if (!colCovered[j] && costMatrix[i][j] < minCost) {
                    minCost = costMatrix[i][j];
                }
            }
        }
    }
    return minCost;
}

bool ckmin(double& a, const double& b) { return b < a ? a = b, 1 : 0; }

// Function to perform the Hungarian algorithm
std::vector<std::vector<int>> hungarianAlgorithm(const std::vector<std::vector<double>>& costMatrix, int k) {
    const int W = costMatrix.size();
    const int H = costMatrix[0].size();
    assert(W <= H);

    // wells[h] = wells assigned to h-th house, or -1 if no wells assigned
    // note: a H-th house was added for convenience
    vector<int> wells(H + 1, -1);
    vector<double> wellPotentials(W), housePotentials(H + 1);  // potentials
    // -housePotentials[H] will equal the sum of all deltas
    vector<vector<int>> answers(W/k);
    const double inf = numeric_limits<double>::max();

    for (int curWell = 0; curWell < W; ++curWell) { // assign curWell-th wells
        int w_cur = H;
        wells[w_cur] = curWell;
        vector<double> min_to(H + 1, inf); // min reduced cost over edges from Z to well h
        vector<int> prv(H + 1, -1);    // previous well on alternating path
        vector<bool> in_Z(H + 1);      // whether well is in Z

        while (wells[w_cur] != -1) {   // runs at most curWell + 1 times
            in_Z[w_cur] = true;
            const int j = wells[w_cur];
            double delta = inf;
            int w_next;
            for (int w = 0; w < H; ++w) {
                if (!in_Z[w]) {
                    if (ckmin(min_to[w], costMatrix[j][w] - wellPotentials[j] - housePotentials[w]))
                        prv[w] = w_cur;
                    if (ckmin(delta, min_to[w])) w_next = w;
                }
            }
            // delta will always be non-negative,
            // except possibly during the first time this loop runs
            // if any entries of C[curWell] are negative
            for (int w = 0; w <= H; ++w) {
                if (in_Z[w]) wellPotentials[wells[w]] += delta, housePotentials[w] -= delta;
                else min_to[w] -= delta;
            }
            w_cur = w_next;
        }
        // update assignments along alternating path
        for (int w; w_cur != H; w_cur = w) {
            wells[w_cur] = wells[w = prv[w_cur]]; 
        }

        if (curWell == W - 1) { // Add assignments only in the last iteration
            for (int h = 0; h < H; ++h) {
                //cout << wells[h] << " ";
                if (wells[h] != -1) {
                    int well_index = wells[h] / k;
                    answers[well_index].push_back(h);
                    // cout << "Well: " << well_index << " wells: " << wells[h] << "\n";
                }
            }
        }
    }

    return answers;
}

// Function to implement the algorithm to find the cheapest water supply
vector<vector<int>> findCheapestWaterSupply(vector<Point> wells, vector<Point> houses, int k) {
    int kn = houses.size(); // Number of houses

    std::vector<Point> duplicatedWells;
    for (Point element : wells) {
        for (int i = 0; i < k; i++) {
            duplicatedWells.push_back(element);
        }
    }

    // Create a matrix to store the costs
    vector<vector<double>> costMatrix(kn, vector<double>(kn, 0));

    // Calculate costs between each well and each wells
    for (int i = 0; i < kn; ++i) {
        for (int j = 0; j < kn; ++j) {
            double d = euclideanDistance(duplicatedWells[i], houses[j]);
            costMatrix[i][j] = d;
        }
    }

    // Convert assignment matrix to vector<vector<int>>
    vector<vector<int>> result = hungarianAlgorithm(costMatrix, k);
    // Hungarian algorithm

    return result;
}


// Custom comparison function for std::next_permutation
bool comparePoints(const Point& p1, const Point& p2) {
    // Comparison based on x-coordinate and then y-coordinate
    return (p1.x < p2.x) || (p1.x == p2.x && p1.y < p2.y);
}

// Function to check the solution
bool checkSolution(vector<vector<int>>& solution, vector<Point>& wells, vector<Point>& houses, int k) {
    double totalCost = 0.0;

    // Calculate total cost of the solution
    for (int i = 0; i < solution.size(); ++i) {
        Point well = wells[i];
        const vector<int>& assignedHouses = solution[i];

        for (int houseIndex : assignedHouses) {
            Point house = houses[houseIndex];
            totalCost += euclideanDistance(well, house);
        }
    }

    // Calculate total cost of all possible valid assignments
    double minCost = numeric_limits<double>::max();
    vector<Point> bestAssignment; // Store houses assignment for the best cost
    sort(houses.begin(), houses.end(), comparePoints); // Sort houses to generate permutations
    do {
        double currentCost = 0.0;
        for (int i = 0; i < wells.size(); ++i) {
            for (int j = 0; j < k; ++j) {
                currentCost += euclideanDistance(wells[i], houses[i * k + j]);
            }
        }
        if (currentCost < minCost) {
            minCost = currentCost;
            bestAssignment = houses;
        }
    } while (next_permutation(houses.begin(), houses.end(), comparePoints));

    const double tolerance = 1e-3;

    // Compare the total cost of the solution with the minimum cost
    cout << "Correct assignments:" << endl;
    for (int i = 0; i < wells.size(); ++i) {
        cout << "Well " << i << " " << "(" << wells[i].x << "," << wells[i].y << ") " << " assigned to houses: ";
        for (int j = 0; j < k; ++j) {
            int houseIndex = i * k + j;
            cout << "(" << bestAssignment[houseIndex].x << "," << bestAssignment[houseIndex].y << ") ";
        }
        cout << endl;
    }
    cout << "Checked assignments:" << endl;
    for (int i = 0; i < solution.size(); ++i) {
        cout << "Well " << i << " " << "(" << wells[i].x << "," << wells[i].y << ") " << " assigned to houses: ";
        for (int houseIndex : solution[i]) {
            Point house = houses[houseIndex];
            cout << "(" << house.x << "," << house.y << ") ";
        }
        cout << endl;
    }
    cout << "Correct total cost: " << minCost << " ";
    cout << "Checked total cost: " << totalCost << endl;

    if ((totalCost - minCost) <= tolerance && solution[0].size() > 0) {
        return true; // Solution is valid
    }
    else {
        return false; // Solution is invalid
    }
}


int readInputFromFile(const string& filename, int& n, int& k, vector<Point>& wells, vector<Point>& houses) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return -1;
    }

    file >> n >> k;
    wells.resize(n);
    houses.resize(n * k);

    for (int i = 0; i < n; ++i) {
        double x, y;
        file >> x >> y;
        wells[i] = { x, y };
    }

    for (int i = 0; i < n * k; ++i) {
        double x, y;
        file >> x >> y;
        houses[i] = { x, y };
    }

    file.close();
    return 0;
}

int main() {
    
    // std::ofstream devnull("/dev/null");
    // std::streambuf* oldbuf = std::cout.rdbuf(devnull.rdbuf());

    string inputDir = "./inputs/"; // Project directory
    if (!fs::exists(inputDir)) {
        fs::create_directory(inputDir);
    }

    int inputFileCount = 0;

    // Perform 1000 tests
    for (int i = 1; i <= 10; ++i) {
        string filename = inputDir + "random_input_" + to_string(i) + ".txt";
        generateRandomInput(filename);

        // Read input from the generated file
        int n, k;
        vector<Point> wells, houses;
        if (readInputFromFile(filename, n, k, wells, houses) == -1) {
            cerr << "Error reading input from file: " << filename << endl;
            continue;
        }

        cout << "File " << filename << ":" << endl;
        cout << "Number of wells: " << n << endl;
        cout << "Number of houses each well supplies: " << k << endl;

        cout << "Wells: " << endl;
        for (const auto& well : wells) {
            cout << well.x << ", " << well.y << endl;
        }

        cout << "Houses: " << endl;
        for (const auto& house : houses) {
            cout << house.x << ", " << house.y << endl;
        }

        // Call the algorithm function
        vector<vector<int>> solution = findCheapestWaterSupply(wells, houses, k);
        
        // unblock
        // std::cout.rdbuf(oldbuf);
        // devnull.close();

        // Check the solution
        if (checkSolution(solution, wells, houses, k)) {
            cout << "Solution " << filename << " is valid." << endl;
        }
        else {
            cout << "Solution " << filename << " is invalid." << endl;
            return -1;
        }

        // block
        // std::ofstream devnull("/dev/null");
        // std::streambuf* oldbuf = std::cout.rdbuf(devnull.rdbuf());


        inputFileCount++;
    }

    if (inputFileCount == 0) {
        cout << "No input files found in directory." << endl;
    }

    return 0;

}
