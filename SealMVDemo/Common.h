#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "seal/seal.h"
using namespace std;
using namespace seal;

void MVMultiplicationSparse(int dimension, int diagonalStart, int diagonalEnd, bool rescale);
void Sequential(int dimension,bool rescale);
void MVMultiplicationDense(int dimension, bool rescale);

// Helper function that prints a matrix (vector of vectors)
template <typename T>
inline void print_full_matrix(vector<vector<T>> matrix, int precision = 3)
{
    // save formatting for cout
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(precision);
    int row_size = matrix.size();
    int col_size = matrix[0].size();
    for (unsigned int i = 0; i < row_size; i++)
    {
        cout << "[";
        for (unsigned int j = 0; j < col_size - 1; j++)
        {
            cout << matrix[i][j] << ", ";
        }
        cout << matrix[i][col_size - 1];
        cout << "]" << endl;
    }
    cout << endl;
    // restore old cout formatting
    cout.copyfmt(old_fmt);
}

// Helper function that prints parts of a matrix (only squared matrix)
template <typename T>
inline void print_partial_matrix(vector<vector<T>> matrix, int print_size = 3, int precision = 3)
{
    // save formatting for cout
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(precision);

    int row_size = matrix.size();
    int col_size = matrix[0].size();

    // Boundary check
    if(row_size < 2 * print_size && col_size < 2 * print_size)
    {
        cerr << "Cannot print matrix with these dimensions: " << to_string(row_size) << "x" << to_string(col_size) << ". Increase the print size" << endl;
        return;
    }
    // print first 4 elements
    for (unsigned int row = 0; row < print_size; row++)
    {
        cout << "\t[";
        for (unsigned int col = 0; col < print_size; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << "..., ";
        for (unsigned int col = col_size - print_size; col < col_size - 1; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << matrix[row][col_size - 1];
        cout << "]" << endl;
    }
    cout << "\t..." << endl;

    for (unsigned int row = row_size - print_size; row < row_size; row++)
    {
        cout << "\t[";
        for (unsigned int col = 0; col < print_size; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << "..., ";
        for (unsigned int col = col_size - print_size; col < col_size - 1; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << matrix[row][col_size - 1];
        cout << "]" << endl;
    }

    cout << endl;
    // restore old cout formatting
    cout.copyfmt(old_fmt);
}

template <typename T>
inline void print_partial_vector(vector<T> vec, int size, int print_size = 3, int precision = 3)
{
    // save formatting for cout
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(precision);

    int row_size = size;

    // Boundary check
    if (row_size < 2 * print_size)
    {
        cerr << "Cannot print vector with these dimensions: " << to_string(row_size) << ". Increase the print size" << endl;
        return;
    }

    cout << "\t[";
    for (unsigned int row = 0; row < print_size; row++)
    {
        cout << vec[row] << ", ";
    }
    cout << "..., ";

    for (unsigned int row = row_size - print_size; row < row_size - 1; row++)
    {
        cout << vec[row] << ", ";
    }
    cout << vec[row_size - 1] << "]\n";

    cout << endl;
    // restore old cout formatting
    cout.copyfmt(old_fmt);
}

// Gets a diagonal from a matrix U
template <typename T>
inline vector<T> get_diagonal(int position, vector<vector<T>> U)
{

    vector<T> diagonal(U.size());

    int k = 0;
    // U(0,l) , U(1,l+1), ... ,  U(n-l-1, n-1)
    for (int i = 0, j = position; (i < U.size() - position) && (j < U.size()); i++, j++)
    {
        diagonal[k] = U[i][j];
        k++;
    }
    for (int i = U.size() - position, j = 0; (i < U.size()) && (j < position); i++, j++)
    {
        diagonal[k] = U[i][j];
        k++;
    }

    return diagonal;
}

inline Ciphertext Linear_Transform_Plain(Ciphertext ct, vector<Plaintext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params, bool rescale)
{
    cout << "    + Invoked Linear_Transform_Plain to compute C_vec . P_mat" << endl;

    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);

    // Fill ct with duplicate
    Ciphertext ct_rot;

    evaluator.rotate_vector(ct, -U_diagonals.size(), gal_keys, ct_rot);
    Ciphertext ct_new;
    evaluator.add(ct, ct_rot, ct_new);
    cout << "    + Scale Linear_Transform_Plain add: " << log2(ct_new.scale()) << " bits" << endl;

    vector<Ciphertext> ct_result(U_diagonals.size());
    evaluator.multiply_plain(ct_new, U_diagonals[0], ct_result[0]);

    for (int l = 1; l < U_diagonals.size(); l++)
    {
        Ciphertext temp_rot;
        evaluator.rotate_vector(ct_new, l, gal_keys, temp_rot);
        evaluator.multiply_plain(temp_rot, U_diagonals[l], ct_result[l]);
    }

    Ciphertext ct_prime;
    evaluator.add_many(ct_result, ct_prime);
    if (rescale) 
    {
        cout << "    + Scale Linear_Transform_Plain add_many: " << log2(ct_prime.scale()) << " bits" << endl;
        evaluator.rescale_to_next_inplace(ct_prime);
        cout << "    + after rescale Linear_Transform_Plain add_many: " << log2(ct_prime.scale()) << " bits" << endl;
    }    

    return ct_prime;
}

inline Ciphertext Linear_Transform_Cipher(Ciphertext ct, vector<Ciphertext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params, bool rescale)
{
    cout << "    + Invoked Linear_Transform_Cipher C_vec . C_mat:" << endl;
    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);

    // Fill ct with duplicate
    Ciphertext ct_rot;

    evaluator.rotate_vector(ct, -U_diagonals.size(), gal_keys, ct_rot);    
    Ciphertext ct_new;
    evaluator.add(ct, ct_rot, ct_new);

    cout << "    + Scale add: " << log2(ct_new.scale()) << " bits" << endl;

    vector<Ciphertext> ct_result(U_diagonals.size());
    evaluator.multiply(ct_new, U_diagonals[0], ct_result[0]);       
    for (int l = 1; l < U_diagonals.size(); l++)
    {
        Ciphertext temp_rot;
        evaluator.rotate_vector(ct_new, l, gal_keys, temp_rot);
        evaluator.multiply(temp_rot, U_diagonals[l], ct_result[l]);              
    } 
    
    Ciphertext ct_prime;
    evaluator.add_many(ct_result, ct_prime);
    if (rescale && log2(ct_prime.scale()) > 40)
    {
        cout << "    + Scale Linear_Transform_Cipher addmany: " << log2(ct_prime.scale()) << " bits" << endl;        
        evaluator.rescale_to_next_inplace(ct_prime);
        cout << "    + after rescale Linear_Transform_Cipher addmany: " << log2(ct_prime.scale()) << " bits" << endl;
    }    

    return ct_prime;
}

inline void test_Linear_Transformation(int dimension, vector<vector<double>> input_matrix, vector<double> input_vec)
{
    vector<double> result(dimension);
    int k = 0;
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            result[k] += input_matrix[i][j] * input_vec[j];
        }
        k++;
    }

    // Print Result vector
    print_partial_vector(result, dimension);
}

//return result
inline vector<double> get_Linear_Transformation_expected_vector(int dimension, vector<vector<double>> input_matrix, vector<double> input_vec)
{
    vector<double> result(dimension);
    int k = 0;
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            result[k] += input_matrix[i][j] * input_vec[j];
        }
        k++;
    }

    return result;
}

//return result
inline vector<double> get_Linear_Transformation_expected_sequentialVector(int dimension, vector<vector<double>> input_matrix)
{
    vector<double> result(dimension);
    vector<double> startInput = input_matrix[0];
    vector<vector<double>> outputvector(dimension);
    int k = 0;
    for (int i = 0; i < dimension; i++)
    {
        if (i == 0)
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, startInput);
        }
        else
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, outputvector[i-1]);
        }               
    }
    int size = outputvector.size();
    return outputvector[size-1];
}

inline void print_error_difference(vector<double> actualResult, vector<double> expectedResult, int dimention)
{    
    for (unsigned int i = 0; i < dimention; i++)
    {
        cout << "index " << i << ": " << expectedResult[i] - actualResult[i] << "\n";
    }
}

inline void get_max_error_norm(vector<double> actualResult, vector<double> expectedResult, int dimension)
{    
    vector<double> error(dimension);
    for (unsigned int i = 0; i < dimension; i++)
    {
        double difference = expectedResult[i] - actualResult[i];
        error[i] = abs(difference);        
    }

    // Find the maximum element in error vector
    cout << "\nMax Error Norm = " << *max_element(error.begin(), error.end()) << "\n";    
    
}
