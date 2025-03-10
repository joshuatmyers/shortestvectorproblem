#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>

// Returns the dot product of two vectors
double dotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2);

// DFS to calculate the size of all possible linear combinations in a given range 
void findShortestVector(const std::vector<std::vector<double>>& basis, 
                        std::vector<int>& coefficients, 
                        size_t depth, 
                        double radius, 
                        std::vector<double>& shortestVector, 
                        double& shortestMagnitude, 
                        std::vector<double>& currentCombination, 
                        double currentMagnitude);

// Orthogonalise the provided basis using Gram-Schmidt
void gramSchmidt(const std::vector<std::vector<double>>& basis, std::vector<std::vector<double>>& mu, std::vector<double>& squaredNorms);

// Check Lovasz condition
bool checkLovaszCondition(const std::vector<double>& squaredNorms, const std::vector<std::vector<double>>& mu, size_t k);

// LLL reduction algorithm
std::vector<std::vector<double>> LLLreduction(std::vector<std::vector<double>>& basis);

// Calculate dynamic radius
void dynamicRadius(const std::vector<std::vector<double>>& basis, double& radius);

#endif