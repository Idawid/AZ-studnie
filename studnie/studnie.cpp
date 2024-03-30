#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <filesystem>

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

// Function to implement the algorithm to find the cheapest water supply
// Function to implement the algorithm to find the cheapest water supply
vector<vector<int>> findCheapestWaterSupply(vector<Point> wells, vector<Point> houses, int k) {
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
        const vector<int>& distancesA = a.second; // well indices
        const vector<int>& distancesB = b.second; // well indices

        // Calculate the difference in distances for houses a and b
        double diffA, diffB;
        if (wells.size() == 1) {
            // If there's only one well, set the difference to 0 for both a and b
            diffA = diffB = 0.0;
        }
        else {
            // Calculate the difference in distances for houses a and b using the two closest wells
            diffA = abs(euclideanDistance(houses[a.first], wells[distancesA[1]]) - euclideanDistance(houses[a.first], wells[distancesA[0]]));
            diffB = abs(euclideanDistance(houses[b.first], wells[distancesB[1]]) - euclideanDistance(houses[b.first], wells[distancesB[0]]));
        }

        // Sort in descending order of the difference in distances
        return diffA > diffB;
    });

    // 3. Greedy Assignment with Prioritization
    vector<vector<int>> solution(wells.size());
    for (const auto& pair : houseWellAssignment) {
        int houseIndex = pair.first;
        const vector<int>& wellIndices = pair.second;

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
}


// Custom comparison function for std::next_permutation
bool comparePoints(const Point& p1, const Point& p2) {
    // Comparison based on x-coordinate and then y-coordinate
    return (p1.x < p2.x) || (p1.x == p2.x && p1.y < p2.y);
}

// Function to check the solution
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

    // Perform 100 tests
    for (int i = 1; i <= 100; ++i) {
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
            cout << "Solution is valid." << endl;
        }
        else {
            cout << "Solution is invalid." << endl;
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
    string directory = "./";
    int inputFileCount = 0;

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            string filename = entry.path().string();
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

            // Check the solution
            if (checkSolution(solution, wells, houses, k)) {
                cout << "Solution is valid." << endl;
            }
            else {
                cout << "Solution is invalid." << endl;
            }

            inputFileCount++;
        }
    }

    if (inputFileCount == 0) {
        cout << "No input files found in directory." << endl;
    }

    return 0;
    */
}
