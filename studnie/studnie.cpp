#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/matching.h>
#include <opencv2/opencv.hpp>

#define E 10000

const int INF = 1e9;

using namespace std;
using namespace lemon;
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

// Function to perform the Hungarian algorithm
vector<vector<int>> hungarianAlgorithm(const vector<vector<double>>& costMatrix, int k) {
    int numRows = costMatrix.size();
    int numCols = costMatrix[0].size();

    // Step 1: Subtract the minimum value of each row from its elements and cover the row
    vector<bool> rowCovered(numRows, false);
    vector<bool> colCovered(numCols, false);
    vector<vector<double>> reducedCostMatrix(costMatrix);

    for (int i = 0; i < numRows; ++i) {
        double minVal = *min_element(reducedCostMatrix[i].begin(), reducedCostMatrix[i].end());
        for (int j = 0; j < numCols; ++j) {
            reducedCostMatrix[i][j] -= minVal;
        }
    }

    // Step 2: Subtract the minimum value of each column from its elements and cover the column
    for (int j = 0; j < numCols; ++j) {
        double minVal = numeric_limits<double>::infinity();
        for (int i = 0; i < numRows; ++i) {
            minVal = min(minVal, reducedCostMatrix[i][j]);
        }
        for (int i = 0; i < numRows; ++i) {
            reducedCostMatrix[i][j] -= minVal;
        }
    }

    // Step 3: Cover all zeros in the matrix using the minimum number of lines
    while (true) {
        // Mark all zeros with a star if no other stars are in its row or column
        vector<vector<bool>> starredZeros(numRows, vector<bool>(numCols, false));
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (reducedCostMatrix[i][j] == 0 && !rowCovered[i] && !colCovered[j]) {
                    starredZeros[i][j] = true;
                    rowCovered[i] = true;
                    colCovered[j] = true;
                    break;
                }
            }
        }

        // Count the number of covered columns
        int numCoveredCols = 0;
        for (bool col : colCovered) {
            if (col) {
                numCoveredCols++;
            }
        }

        // If the number of covered columns equals the number of rows, we're done
        if (numCoveredCols == numRows) {
            break;
        }

        // Find an uncovered zero and prime it
        int row, col;
        bool zeroFound = false;
        for (int i = 0; i < numRows && !zeroFound; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (reducedCostMatrix[i][j] == 0 && !rowCovered[i] && !colCovered[j]) {
                    row = i;
                    col = j;
                    zeroFound = true;
                    break;
                }
            }
        }

        // If no uncovered zero was found, proceed to Step 4
        if (!zeroFound) {
            // Find the minimum uncovered element
            double minUncovered = findMinCost(reducedCostMatrix, rowCovered, colCovered);

            // Subtract the minimum uncovered element from all uncovered rows, and add it to all covered columns
            for (int i = 0; i < numRows; ++i) {
                if (rowCovered[i]) {
                    for (int j = 0; j < numCols; ++j) {
                        if (colCovered[j]) {
                            reducedCostMatrix[i][j] += minUncovered;
                        }
                    }
                }
                else {
                    for (int j = 0; j < numCols; ++j) {
                        if (!colCovered[j]) {
                            reducedCostMatrix[i][j] -= minUncovered;
                        }
                    }
                }
            }
        }
        else {
            // Prime the found zero
            starredZeros[row][col] = true;

            // Find a starred zero in the same row
            bool foundStarredZero = false;
            for (int j = 0; j < numCols && !foundStarredZero; ++j) {
                if (starredZeros[row][j] && j != col) {
                    col = j;
                    foundStarredZero = true;
                }
            }

            // If a starred zero was found in the same row, cover the row and uncover the column containing the starred zero
            if (foundStarredZero) {
                rowCovered[row] = true;
                colCovered[col] = false;
            }
            else {
                // Otherwise, proceed to Step 3 with the new zero
                colCovered[col] = true;
            }
        }
    }

    // Construct the assignment vector
    vector<vector<int>> assignment(costMatrix.size() / k);
    for (int j = 0; j < numCols; ++j) {
        for (int i = 0; i < numRows; ++i) {
            if (reducedCostMatrix[i][j] == 0 && !colCovered[j]) {
                assignment[j / k].push_back(i / k); // Map back to original wells
                colCovered[j] = true;
                break;
            }
        }
    }

    return assignment;
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

    // Calculate costs between each well and each house
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

    /*
    // 1. House Perspective
    vector<pair<int, vector<int>>> houseWellAssignment; // <house index, well indices sorted from closest>
    for (int i = 0; i < houses.size(); ++i) {
        Point house = houses[i];
        vector<pair<int, double>> distances; // Pair of (well index, distance to the well)

        // Calculate distance to all wells
        for (int j = 0; j < wells.size(); ++j) {
            double distance = euclideanDistance(house, wells[j]);
            distances.push_back({ j, distance });
        }

        // Sort distances in ascending order
        sort(distances.begin(), distances.end(), [](const auto& a, const auto& b) {
            return a.second < b.second;
            });

        // Store house index and sorted well indices
        vector<int> wellIndices;
        for (const auto& pair : distances) {
            wellIndices.push_back(pair.first);
        }
        houseWellAssignment.push_back({ i, wellIndices });
    }

    // 2. Prioritization
    sort(houseWellAssignment.begin(), houseWellAssignment.end(), [&](const auto& a, const auto& b) {
        const vector<int>& wellsA = a.second; // well indices
        const vector<int>& wellsB = b.second; // well indices

        // Calculate the difference in distances for houses a and b
        double diffA, diffB;
        if (wells.size() == 1) {
            // If there's only one well, set the difference to 0 for both a and b
            diffA = diffB = 0.0;
        }
        else {
            // Calculate the difference in distances for houses a and b using the two closest wells
            diffA = abs(euclideanDistance(houses[a.first], wells[wellsA[1]]) - euclideanDistance(houses[a.first], wells[wellsA[0]]));
            diffB = abs(euclideanDistance(houses[b.first], wells[wellsB[1]]) - euclideanDistance(houses[b.first], wells[wellsB[0]]));
        }

        // Sort in descending order of the difference in distances
        return diffA > diffB;
    });

    // 3. Greedy Assignment with Prioritization
    vector<vector<int>> solution(wells.size());

    /*for (int i = 0; i < houses.size(); i++) {
        if (!houseWellAssignment.empty()) {
            sort(houseWellAssignment.begin(), houseWellAssignment.end(), [&](const auto& a, const auto& b) {
                const vector<int>& wellsA = a.second; // well indices
                const vector<int>& wellsB = b.second; // well indices

                // Calculate the difference in distances for houses a and b
                double diffA, diffB;
                if (wellsA.size() == 1 || wellsB.size() == 1) {
                    // If there's only one well, set the difference to 0 for both a and b
                    diffA = diffB = 0.0;
                }
                else {
                    // Calculate the difference in distances for houses a and b using the two closest wells
                    diffA = abs(euclideanDistance(houses[a.first], wells[wellsA[1]]) - euclideanDistance(houses[a.first], wells[wellsA[0]]));
                    diffB = abs(euclideanDistance(houses[b.first], wells[wellsB[1]]) - euclideanDistance(houses[b.first], wells[wellsB[0]]));
                }

                // Sort in descending order of the difference in distances
                return diffA > diffB;
                });

            int houseIndex = houseWellAssignment[0].first;
            int closestWellIndex = (houseWellAssignment[0].second)[0];

            cout << "Connecting house " << houseIndex << " (" << houses[houseIndex].x << "," << houses[houseIndex].y << ") " <<
                "with well " << closestWellIndex << " (" << wells[closestWellIndex].x << "," << wells[closestWellIndex].y << ") " "\n";

            solution[closestWellIndex].push_back(houseIndex);

            // Erase the element from houseWellAssignment that contains houseIndex
            houseWellAssignment.erase(std::remove_if(houseWellAssignment.begin(), houseWellAssignment.end(),
                [houseIndex](const std::pair<int, std::vector<int>>& pair) {
                    return pair.first == houseIndex;
                }),
                houseWellAssignment.end());

            // Delete the closest well index for all houses
            if (solution[closestWellIndex].size() == k) {
                for_each(houseWellAssignment.begin(), houseWellAssignment.end(),
                    [closestWellIndex](pair<int, vector<int>>& pair) {
                        pair.second.erase(remove(pair.second.begin(), pair.second.end(), closestWellIndex), pair.second.end());
                    });
            }
        }
    }
    */

    /*
    for (const auto& pair : houseWellAssignment) {
        int houseIndex = pair.first;
        const vector<int>& wellIndices = pair.second;

        sort(houseWellAssignment.begin(), houseWellAssignment.end(), [&](const auto& a, const auto& b) {
            const vector<int>& wellsA = a.second; // well indices
            const vector<int>& wellsB = b.second; // well indices

            // Calculate the difference in distances for houses a and b
            double diffA, diffB;
            if (wells.size() == 1) {
                // If there's only one well, set the difference to 0 for both a and b
                diffA = diffB = 0.0;
            }
            else {
                // Calculate the difference in distances for houses a and b using the two closest wells
                diffA = abs(euclideanDistance(houses[a.first], wells[wellsA[1]]) - euclideanDistance(houses[a.first], wells[wellsA[0]]));
                diffB = abs(euclideanDistance(houses[b.first], wells[wellsB[1]]) - euclideanDistance(houses[b.first], wells[wellsB[0]]));
            }

            // Sort in descending order of the difference in distances
            return diffA > diffB;
        });

        // Iterate through sorted well indices
        for (int wellIndex : wellIndices) {
            // Check if the well has reached its capacity
            if (solution[wellIndex].size() < k) {
                // Assign the house index to this well
                solution[wellIndex].push_back(houseIndex);
                break;
            }
        }
    }
    
    // 4. Cost Calculation
    double totalCost = 0.0;

    for (int i = 0; i < solution.size(); ++i) {
        Point well = wells[i];
        const vector<int>& assignedHouses = solution[i];

        cout << "Costs from well " << i << " to assigned houses:" << endl;

        for (int houseIndex : assignedHouses) {
            Point house = houses[houseIndex];
            double distance = euclideanDistance(well, house);
            totalCost += distance;
            cout << "   House " << houseIndex << ": " << distance << endl;
        }
    }

    // Print total cost
    cout << "Total cost of the assignment: " << totalCost << endl;

    return solution;
    */
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
    if ((totalCost - minCost) <= tolerance) {
        return true; // Solution is valid
    }
    else {
        cout << "Solution is invalid." << endl;
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
        cout << "Correct total cost: " << minCost << endl;
        cout << "Checked total cost: " << totalCost << endl;
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
    
    std::ofstream devnull("/dev/null");
    std::streambuf* oldbuf = std::cout.rdbuf(devnull.rdbuf());

    string directory = "./"; // Project directory
    int inputFileCount = 0;

    // Perform 1000 tests
    for (int i = 1; i <= 1000; ++i) {
        string filename = "random_input_" + to_string(i) + ".txt";
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
        std::cout.rdbuf(oldbuf);
        devnull.close();

        // Check the solution
        if (checkSolution(solution, wells, houses, k)) {
            cout << "Solution " << filename << " is valid." << endl;
        }
        else {
            cout << "Solution " << filename << " is invalid." << endl;
            return -1;
        }

        // block
        std::ofstream devnull("/dev/null");
        std::streambuf* oldbuf = std::cout.rdbuf(devnull.rdbuf());


        inputFileCount++;
    }

    if (inputFileCount == 0) {
        cout << "No input files found in directory." << endl;
    }

    return 0;

    
    /*
    string filename = "input.txt";

    int n, k;
    vector<Point> wells, houses;

    if (readInputFromFile(filename, n, k, wells, houses) == -1) {
        cerr << "Error reading input from file: " << filename << endl;
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

    // Check the solution
    if (checkSolution(solution, wells, houses, k)) {
        cout << "Solution is valid." << endl;
    }
    else {
        cout << "Solution is invalid." << endl;
    }
    */
    return 0;
}
