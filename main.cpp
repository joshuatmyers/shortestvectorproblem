#include "main.h"

// Returns the dot product of two vectors
double dotProduct(const std::vector<double>& vector1,
                  const std::vector<double>& vector2) {
    if (vector1.size() != vector2.size()) {
        return 0;
    }

    double res = 0;

    for (size_t i = 0; i < vector1.size(); i++) {
        res += vector1[i] * vector2[i];
    }

    return res;
}

// DFS to calculate the size of all linear combinations in a given range
void findShortestVector(const std::vector<std::vector<double>>& basis,
                        std::vector<int>& coefficients,
                        size_t depth,
                        double radius,
                        std::vector<double>& shortestVector,
                        double& shortestMagnitude,
                        std::vector<double>& currentCombination,
                        double currentMagnitude) {

    /*
    basis: the reduced basis that are being used to calculate the shortest vector
    coefficients: store the coefficients of the linear combination
    depth: the number of vectors in the provided basis currently being considered, used for base case
    radius: defines maximum radius to search for SVP, the ceiling of this value if used for the max coefficient
    currentCombination: the result of the linear combination, used for storing the shortest vector
    */

    // If the current magnitude already exceeds the radius, terminate
    if (currentMagnitude > radius) {
        return;
    }

    // Base case when all basis vectors are being considered
    if (depth == basis.size()) {
        // only update when a smaller magnitude is found, which is non-zero
        if (currentMagnitude < shortestMagnitude && currentMagnitude > 0) {
            shortestMagnitude = currentMagnitude;
            shortestVector = currentCombination;
        }
        return;
    }

    // Loop through all possible coefficients iteratively with DFS
    // Handle case that integer is not used for radius
    int maxCoefficient = std::ceil(radius);
    for (int i = -maxCoefficient; i <= maxCoefficient; ++i) {
        coefficients[depth] = i;

        // Update the current combination and its magnitude
        std::vector<double> newCombination = currentCombination;
        double newMagnitude = 0.0;

        // Update new combination with the new coefficient, update magnitude
        for (size_t j = 0; j < basis[depth].size(); ++j) {
            newCombination[j] += i * basis[depth][j];
            newMagnitude += newCombination[j] * newCombination[j];
        }
        newMagnitude = std::sqrt(newMagnitude);

        // Recursively call this function for each combination at current depth
        findShortestVector(basis, coefficients, depth + 1, radius,
                           shortestVector, shortestMagnitude,
                           newCombination, newMagnitude);
    }
}

// Orthogonalise the provided basis using Gram-Schmidt
void gramSchmidt(const std::vector<std::vector<double>>& basis,
                 std::vector<std::vector<double>>& mu,
                 std::vector<double>& squaredNorms) {
    /*
    basis: the basis being orthogonalised
    mu: the list of calculated projections passed to the function to be updated
    squaredNorms: stored euclidean norms of orthogonal basis to be used for Lovascz condition
    */
    size_t n = basis.size();
    std::vector<std::vector<double>> orthogonal(n, std::vector<double>
                                                (basis[0].size(), 0.0));
    // First vector is just normalized
    orthogonal[0] = basis[0];
    // Calculate initial squared euclidean norm
    squaredNorms[0] = dotProduct(orthogonal[0], orthogonal[0]);

    for (size_t i = 1; i < n; ++i) {
        orthogonal[i] = basis[i];
        for (size_t j = 0; j < i; ++j) {
            // Subtract projections to get orthogonal basis
            mu[i][j] = dotProduct(basis[i], orthogonal[j]) / squaredNorms[j];
            for (size_t k = 0; k < orthogonal[i].size(); ++k) {
                orthogonal[i][k] -= mu[i][j] * orthogonal[j][k];
            }
        }
        // Update squared norms; used later to check for Lovascz condition
        squaredNorms[i] = dotProduct(orthogonal[i], orthogonal[i]);
    }
}


// Lovascz condition adapted from https://github.com/orisano/olll/blob/master/olll.py
bool checkLovaszCondition(std::vector<double>& squaredNorms,
                          std::vector<std::vector<double>>& mu,
                          size_t k) {
    return squaredNorms[k] >= ((3.0 / 4.0) - mu[k][k-1] * mu[k][k-1]) * squaredNorms[k-1];
}

// LLL algorithm - adapted from https://math.mit.edu/~apost/courses/18.204-2016/18.204_Xinyue_Deng_final_paper.pdf
std::vector<std::vector<double>> LLLreduction(std::vector<std::vector<double>>& basis) {
    size_t n = basis.size();
    // Store projections and squared norms for Lovasz condition
    std::vector<std::vector<double>> mu(n, std::vector<double>(n, 0.0));
    std::vector<double> squaredNorms(n, 0.0);

    // Initial Gram-Schmidt orthogonalisation
    gramSchmidt(basis, mu, squaredNorms);

    size_t k = 1;
    while (k < n) {
        // Size condition and reduction
        for (int j = k - 1; j >= 0; --j) {
            double m = mu[k][j];
            if (m > 0.5) {
                m = std::round(m);
                for (size_t l = 0; l < basis[k].size(); ++l) {
                    basis[k][l] -= m * basis[j][l];
                }
                // Re-orthogonalize after modification
                gramSchmidt(basis, mu, squaredNorms);
            }
        }

        // Step 2: Check Lovasz condition
        if (checkLovaszCondition(squaredNorms, mu, k)) {
            k = k + 1;
        } else {
            // Swap vectors if condition isn't met
            std::swap(basis[k], basis[k - 1]);
            // Re-orthogonalize after swap
            gramSchmidt(basis, mu, squaredNorms);
            // Update value of k, cannot go below 1
            k = std::max<size_t>(k - 1, 1);
        }
    }
    return basis;
}

void dynamicRadius(const std::vector<std::vector<double>>& basis,
                   double& radius) {
    radius = 0;
    for (size_t i = 0; i < basis.size(); i++) {
        double norm = 0;
        for (size_t j = 0; j < basis[0].size(); j++) {
            norm += basis[i][j] * basis[i][j];
        }
        radius += std::sqrt(norm);
    }
    radius /= basis.size();
}

int main(int argc, char* argv[]) {
    // Check if at least one basis vector is provided
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [basis_vector1] [basis_vector2] ... [basis_vectorN]" << std::endl;
        return 1;
    }

    // Vector to store basis vectors
    // define a vector that contains vectors which contain float values
    std::vector<std::vector<double>> basisVectors;
    std::vector<double> basisVector;
    bool vectorEnd = false;
    
    std::string delimiter = " ";

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s.find('[') != std::string::npos) {
            // clear previous vector for the start of a new vector
            basisVector.clear();
            // remove square brackets from string
            s.erase(std::remove(s.begin(), s.end(), '['), s.end());
        }

        if (s.find(']') != std::string::npos) {
            s.erase(std::remove(s.begin(), s.end(), ']'), s.end());
            // update boolean to refer to later
            vectorEnd = true;
        }

        size_t pos = 0;
        std::string token;

        // modified method from Vincenzo Pii:
        // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            double value = std::stof(token);
            basisVector.push_back(value);
            s.erase(0, pos + delimiter.length());
        }

        // Convert the last token to float and add to the vector
        double lastValue = std::stof(s);
        basisVector.push_back(lastValue);

        // Add the basis vector to the vector of basis vectors
        // check for boolean if this argument contains a closing square bracket
        if (vectorEnd == true) {
            // std::cout << "Vector Pushed" << std::endl;
            basisVectors.push_back(basisVector);
            vectorEnd = false;
        }
    }

    // Start measuring time for reduction and enumeration
    auto start = std::chrono::high_resolution_clock::now();


    std::vector<std::vector<double>> lovasczBasis = LLLreduction(basisVectors);

    std::cout << "Lovascz Reduced Basis: "<< std::endl;
    for (size_t i = 0; i < lovasczBasis.size(); i++) {
        for (size_t j = 0; j < lovasczBasis[i].size(); j++)
            std::cout << lovasczBasis[i][j] << " ";
        std::cout << std::endl;
    }

    // Enumeration of the lovascz basis
    // create dynamic radius size for different basis
    double radius;

    dynamicRadius(lovasczBasis, radius);

    // store coefficients for linar combinations
    std::vector<int> coefficients(lovasczBasis.size(), 0);
    std::vector<double> shortestVector;
    // initialise shortestMagnitude to maximum possible value for a double
    double shortestMagnitude = std::numeric_limits<double>::max();

    std::cout << "Finding the shortest vector... " << std::endl;

    std::vector<double> currentCombination(lovasczBasis[0].size(), 0.0);
    double currentMagnitude = 0.0;
    findShortestVector(lovasczBasis, coefficients, 0, radius,
                       shortestVector, shortestMagnitude,
                       currentCombination, currentMagnitude);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast
                    <std::chrono::milliseconds>(stop - start);

    // Output the shortest vector
    if (!shortestVector.empty()) {
        std::cout << "The shortest vector within radius " << radius << " is: ";
        for (size_t i = 0; i < shortestVector.size(); i++) {
            std::cout << shortestVector[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "The size of the shortest vector: " <<
        std::sqrt(dotProduct(shortestVector, shortestVector)) << std::endl;
        std::cout << "Time taken for calculations : " <<
        duration.count() << " milliseconds" << std::endl;

        // Write results to .txt file
        std::ofstream outputFile("result.txt");
        if (outputFile.is_open()) {
            // write data to the file
            outputFile << std::sqrt(dotProduct(shortestVector, shortestVector));
            outputFile.close(); // close the file when done
            std::cout << "Data was written to output.txt\n";
        }
        else {
            std::cerr << "Error opening file\n";
        }
    } else {
        std::cout << "No vector found within the specified radius.\n" 
        << std::endl;
    }

    return 0;
}