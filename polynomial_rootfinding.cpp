#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <sstream>

// Ensure M_PI is defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double TOLERANCE = 1e-8;
const int MAX_ITERATIONS = 1000;

// Evaluate the polynomial P(z) using Horner’s method.
std::complex<double> evaluatePolynomial(const std::vector<double> &coeffs, std::complex<double> z)
{
    std::complex<double> result = 0;
    for (double coeff : coeffs)
    {
        result = result * z + coeff;
    }
    return result;
}

// Synthetic division: divides poly by (x - r) and returns quotient in 'quotient'.
// Returns true if the remainder is within tolerance (i.e. r is a root).
bool syntheticDivision(const std::vector<double> &poly, double r, std::vector<double> &quotient)
{
    int n = poly.size() - 1;
    quotient.resize(n);
    double carry = poly[0];
    quotient[0] = carry;
    for (int i = 1; i < n; i++)
    {
        carry = carry * r + poly[i];
        quotient[i] = carry;
    }
    double remainder = carry * r + poly[n];
    return std::abs(remainder) < TOLERANCE;
}

// Extracts integer roots from the polynomial (if any) using a simple rational root test.
// When a root is found, it is factored out (possibly repeatedly) via synthetic division.
void extractIntegerRoots(std::vector<double> &poly, std::vector<std::complex<double>> &roots)
{
    bool foundRoot = true;
    while (foundRoot && poly.size() > 1)
    {
        foundRoot = false;
        double constantTerm = poly.back();
        // Special case: constant term 0 means x=0 is a root.
        if (std::abs(constantTerm) < TOLERANCE)
        {
            roots.push_back(0);
            std::vector<double> newPoly(poly.size() - 1);
            for (size_t i = 0; i < newPoly.size(); i++)
            {
                newPoly[i] = poly[i];
            }
            poly = newPoly;
            foundRoot = true;
            continue;
        }
        // Check if the constant term is nearly an integer.
        double roundedConst = std::round(constantTerm);
        if (std::abs(constantTerm - roundedConst) < TOLERANCE)
        {
            int intConst = static_cast<int>(std::abs(roundedConst));
            std::vector<int> candidates;
            // Gather all divisors of the constant term.
            for (int i = 1; i <= intConst; i++)
            {
                if (intConst % i == 0)
                {
                    candidates.push_back(i);
                    candidates.push_back(-i);
                }
            }
            // Remove duplicate candidates.
            std::sort(candidates.begin(), candidates.end());
            candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
            // Test each candidate using synthetic division.
            for (int candidate : candidates)
            {
                std::vector<double> quotient;
                if (syntheticDivision(poly, candidate, quotient))
                {
                    // Candidate is a root. Factor it out and add it to the list.
                    roots.push_back(candidate);
                    poly = quotient;
                    foundRoot = true;
                    break; // Restart search on the reduced polynomial.
                }
            }
        }
    }
}

// Durand-Kerner method to find all (complex) roots of the polynomial given by 'coeffs'.
std::vector<std::complex<double>> findPolynomialRoots(const std::vector<double> &coeffs)
{
    int degree = coeffs.size() - 1;
    std::vector<std::complex<double>> roots(degree);
    // Initialize guesses uniformly around a circle (radius can be adjusted).
    double radius = 0.4;
    for (int i = 0; i < degree; ++i)
    {
        double theta = 2 * M_PI * i / degree;
        roots[i] = std::polar(radius, theta);
    }

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter)
    {
        bool converged = true;
        std::vector<std::complex<double>> newRoots = roots;
        for (int i = 0; i < degree; ++i)
        {
            std::complex<double> numerator = evaluatePolynomial(coeffs, roots[i]);
            std::complex<double> denominator = 1.0;
            for (int j = 0; j < degree; ++j)
            {
                if (i != j)
                {
                    denominator *= (roots[i] - roots[j]);
                }
            }
            std::complex<double> delta = numerator / denominator;
            newRoots[i] = roots[i] - delta;
            if (std::abs(delta) > TOLERANCE)
            {
                converged = false;
            }
        }
        roots = newRoots;
        if (converged)
        {
            break;
        }
    }
    return roots;
}

// Merge duplicate (or nearly duplicate) roots that differ by less than TOLERANCE.
std::vector<std::complex<double>> mergeDuplicateRoots(const std::vector<std::complex<double>> &roots)
{
    std::vector<std::complex<double>> merged;
    std::vector<bool> mergedFlag(roots.size(), false);
    for (size_t i = 0; i < roots.size(); ++i)
    {
        if (mergedFlag[i])
            continue;
        std::complex<double> sum = roots[i];
        int count = 1;
        for (size_t j = i + 1; j < roots.size(); ++j)
        {
            if (std::abs(roots[i] - roots[j]) < TOLERANCE)
            {
                sum += roots[j];
                count++;
                mergedFlag[j] = true;
            }
        }
        merged.push_back(sum / static_cast<double>(count));
    }
    return merged;
}

// Function to print the roots.
void printRoots(const std::vector<std::complex<double>> &roots)
{
    std::cout << "Found roots:\n";
    for (const auto &root : roots)
    {
        std::cout << root << "\n";
    }
}

// Interactive input function to read polynomial coefficients.
std::vector<double> getPolynomialFromUser()
{
    int degree;

    while (true)
    {
        std::cout << "\nEnter the degree of the polynomial (must be >= 1): ";
        std::cin >> degree;

        if (std::cin.fail() || degree < 1)
        {
            std::cin.clear();             // Clear the error flag
            std::cin.ignore(10000, '\n'); // Discard invalid input
            std::cout << "Invalid input. Please enter an integer greater than or equal to 1.\n";
        }
        else
        {
            break;
        }
    }

    std::vector<double> poly(degree + 1);

    std::cout << "\nNow enter the " << degree + 1 << " coefficients of the polynomial:\n";
    std::cout << "FORMAT: Enter coefficients in descending order of powers.\n";
    std::cout << "Example: For P(z) = 2z^3 - 4z + 5, enter:\n2 0 -4 5\n\n";

    while (true)
    {
        std::cout << "Enter " << degree + 1 << " numbers (separated by spaces): ";
        std::cin.ignore(); // Clear newline character from previous input
        std::string inputLine;
        std::getline(std::cin, inputLine);
        std::stringstream ss(inputLine);

        bool validInput = true;
        for (int i = 0; i <= degree; i++)
        {
            if (!(ss >> poly[i]))
            {
                validInput = false;
                break;
            }
        }

        if (validInput)
        {
            break;
        }
        else
        {
            std::cout << "Invalid input. Please enter exactly " << degree + 1 << " numeric values.\n";
        }
    }

    return poly;
}

int main()
{
    std::cout << "Polynomial Root Finder\n";
    std::cout << "1. Use a predefined example (P(z) = z^3 - 1)\n";
    std::cout << "2. Enter your own polynomial\n";
    std::cout << "Choose an option: ";

    int option;
    std::cin >> option;

    std::vector<double> coefficients;
    if (option == 1)
    {
        // Predefined polynomial: P(z) = z^3 - 1
        coefficients = {1, 0, 0, -1};
        std::cout << "Using polynomial: z^3 - 1\n";
    }
    else if (option == 2)
    {
        coefficients = getPolynomialFromUser();
    }
    else
    {
        std::cerr << "Invalid option. Exiting.\n";
        return 1;
    }

    std::vector<std::complex<double>> rootsFound;
    // Make a working copy of the coefficients.
    std::vector<double> poly = coefficients;

    // First, extract any obvious integer roots.
    extractIntegerRoots(poly, rootsFound);

    // If there is any polynomial left (degree >= 1), use Durand-Kerner to find its roots.
    if (poly.size() > 1)
    {
        std::vector<std::complex<double>> dkRoots = findPolynomialRoots(poly);
        std::vector<std::complex<double>> mergedDKRoots = mergeDuplicateRoots(dkRoots);
        rootsFound.insert(rootsFound.end(), mergedDKRoots.begin(), mergedDKRoots.end());
    }

    printRoots(rootsFound);

    return 0;
}
