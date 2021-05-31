//#include "seal/seal.h"
//#include<iostream>
//
//using namespace std;
//using namespace seal;
//
//int main()
//{
//	EncryptionParameters parms(scheme_type::ckks);
//	return 0;
//} 

#include <iostream>
#include <iomanip>
#include <fstream>
#include "seal/seal.h"


using namespace std;
using namespace seal;

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
    if (row_size < 2 * print_size && col_size < 2 * print_size)
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
vector<T> get_diagonal(int position, vector<vector<T>> U)
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

Ciphertext Linear_Transform_Plain(Ciphertext ct, vector<Plaintext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params)
{
    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);

    // Fill ct with duplicate
    Ciphertext ct_rot;
    evaluator.rotate_vector(ct, -U_diagonals.size(), gal_keys, ct_rot);
    
    // cout << "U_diagonals.size() = " << U_diagonals.size() << endl;
    Ciphertext ct_new;
    evaluator.add(ct, ct_rot, ct_new);

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

    return ct_prime;
}

Ciphertext Linear_Transform_Cipher(Ciphertext ct, vector<Ciphertext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params)
{
    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);

    // Fill ct with duplicate
    Ciphertext ct_rot;
    evaluator.rotate_vector(ct, -U_diagonals.size(), gal_keys, ct_rot);
    
    // cout << "U_diagonals.size() = " << U_diagonals.size() << endl;
    Ciphertext ct_new;
    evaluator.add(ct, ct_rot, ct_new);

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

    return ct_prime;
}

void test_Linear_Transformation(int dimension, vector<vector<double>> input_matrix, vector<double> input_vec)
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


void TimesForDifferentSets(size_t poly_modulus_degree)
{
    EncryptionParameters params(scheme_type::ckks);
    params.set_poly_modulus_degree(poly_modulus_degree);
    cout << "MAX BIT COUNT: " << CoeffModulus::MaxBitCount(poly_modulus_degree) << endl;
    params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    auto context = SEALContext::SEALContext(params);

    // Generate keys, encryptor, decryptor and evaluator
   /* KeyGenerator keygen(context);
    PublicKey pk = keygen.create_public_key();
    SecretKey sk = keygen.secret_key();
    GaloisKeys gal_keys = keygen.create_galois_keys();*/
    KeyGenerator keygen(context);
    auto sk = keygen.secret_key();
    PublicKey pk;
    keygen.create_public_key(pk);
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    Encryptor encryptor(context, pk);
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);

    // Create CKKS encoder
    CKKSEncoder ckks_encoder(context);
    // Create scale
    cout << "Coeff Modulus Back Value: " << params.coeff_modulus().back().value() << endl;
    double scale = pow(2.0, 40);

    // Set output file
    string filename = "linear_transf_" + to_string(poly_modulus_degree) + ".dat";
    ofstream outf(filename);

    // Handle file error
    if (!outf)
    {
        cerr << "Couldn't open file: " << filename << endl;
        exit(1);
    }

    // Set output script
    string script = "script_linear_transf_" + to_string(poly_modulus_degree) + ".p";
    ofstream outscript(script);

    // Handle script error
    if (!outscript)
    {
        cerr << "Couldn't open file: " << script << endl;
        exit(1);
    }

    // Write to Script
    outscript << "# Set the output terminal" << endl;
    outscript << "set terminal canvas" << endl;
    outscript << "set output \"canvas_linear_transf_" << to_string(poly_modulus_degree) << ".html\"" << endl;
    outscript << "set title \"Linear Transformation Benchmark " << to_string(poly_modulus_degree) << "\"" << endl;
    outscript << "set xlabel 'Dimension'" << endl;
    outscript << "set ylabel 'Time (microseconds)'" << endl;
    outscript << "set logscale" << endl;
    outscript << "set ytics nomirror" << endl;
    outscript << "set xtics nomirror" << endl;
    outscript << "set grid" << endl;
    outscript << "set key outside" << endl;

    outscript << "\n# Set the styling " << endl;
    outscript << "set style line 1\\\n"
        << "linecolor rgb '#0060ad'\\\n"
        << "linetype 1 linewidth 2\\\n"
        << "pointtype 7 pointsize 1.5\n"
        << endl;

    outscript << "set style line 2\\\n"
        << "linecolor rgb '#dd181f'\\\n"
        << "linetype 1 linewidth 2\\\n"
        << "pointtype 5 pointsize 1.5\n"
        << endl;

    outscript << "\nplot 'linear_transf_" << to_string(poly_modulus_degree) << ".dat' index 0 title \"C_Vec * P_Mat\" with linespoints ls 1, \\\n"
        << "'' index 1 title \"C_Vec * C_Mat\"  with linespoints ls 2";
    // Close script
    outscript.close();

    // -------------- FIRST SET ----------------
    int dimension1 = 10;
    cout << "Dimension Set 1: " << dimension1 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set1(dimension1, vector<double>(dimension1));
    vector<vector<double>> pod_matrix2_set1(dimension1, vector<double>(dimension1));

    // Fill input matrices
    double filler = 1.0;
    // Set 1
    double r = ((double)rand() / (RAND_MAX));

    for (int i = 0; i < dimension1; i++)
    {
        for (int j = 0; j < dimension1; j++)
        {
            pod_matrix1_set1[i][j] = r;
            r = ((double)rand() / (RAND_MAX));
            pod_matrix2_set1[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set1);
    print_partial_matrix(pod_matrix2_set1);

    // --------------- SECOND SET ------------------
    int dimension2 = 20;
    cout << "Dimension Set 2: " << dimension2 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set2(dimension2, vector<double>(dimension2));
    vector<vector<double>> pod_matrix2_set2(dimension2, vector<double>(dimension2));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension2; i++)
    {
        for (int j = 0; j < dimension2; j++)
        {
            pod_matrix1_set2[i][j] = filler;
            pod_matrix2_set2[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set2);
    print_partial_matrix(pod_matrix2_set2);

    // --------------- THIRD SET ------------------'
    int dimension3 = 30;
    cout << "Dimension Set 3: " << dimension3 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set3(dimension3, vector<double>(dimension3));
    vector<vector<double>> pod_matrix2_set3(dimension3, vector<double>(dimension3));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension3; i++)
    {
        for (int j = 0; j < dimension3; j++)
        {
            pod_matrix1_set3[i][j] = filler;
            pod_matrix2_set3[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set3);
    print_partial_matrix(pod_matrix2_set3);


    // --------------- FOURTH SET ------------------
    int dimension4 = 40;
    cout << "Dimension Set 4: " << dimension4 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set4(dimension4, vector<double>(dimension4));
    vector<vector<double>> pod_matrix2_set4(dimension4, vector<double>(dimension4));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension4; i++)
    {
        for (int j = 0; j < dimension4; j++)
        {
            pod_matrix1_set4[i][j] = filler;
            pod_matrix2_set4[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set4);
    print_partial_matrix(pod_matrix2_set4);

    // --------------- FIFTH SET ------------------
    int dimension5 = 50;
    cout << "Dimension Set 5: " << dimension5 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set5(dimension5, vector<double>(dimension5));
    vector<vector<double>> pod_matrix2_set5(dimension5, vector<double>(dimension5));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension5; i++)
    {
        for (int j = 0; j < dimension5; j++)
        {
            pod_matrix1_set5[i][j] = filler;
            pod_matrix2_set5[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set5);
    print_partial_matrix(pod_matrix2_set5);

    // --------------- SIXTH SET ------------------
    int dimension6 = 60;
    cout << "Dimension Set 6: " << dimension6 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set6(dimension6, vector<double>(dimension6));
    vector<vector<double>> pod_matrix2_set6(dimension6, vector<double>(dimension6));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension6; i++)
    {
        for (int j = 0; j < dimension6; j++)
        {
            pod_matrix1_set6[i][j] = filler;
            pod_matrix2_set6[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set6);
    print_partial_matrix(pod_matrix2_set6);

    // --------------- SEVENTH SET ------------------
    int dimension7 = 70;
    cout << "Dimension Set 7: " << dimension7 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set7(dimension7, vector<double>(dimension7));
    vector<vector<double>> pod_matrix2_set7(dimension7, vector<double>(dimension7));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension7; i++)
    {
        for (int j = 0; j < dimension7; j++)
        {
            pod_matrix1_set7[i][j] = filler;
            pod_matrix2_set7[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set7);
    print_partial_matrix(pod_matrix2_set7);

    // --------------- EIGHTH SET ------------------
    int dimension8 = 80;
    cout << "Dimension Set 8: " << dimension8 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set8(dimension8, vector<double>(dimension8));
    vector<vector<double>> pod_matrix2_set8(dimension8, vector<double>(dimension8));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension8; i++)
    {
        for (int j = 0; j < dimension8; j++)
        {
            pod_matrix1_set8[i][j] = filler;
            pod_matrix2_set8[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set8);
    print_partial_matrix(pod_matrix2_set8);

    // --------------- NINTH SET ------------------
    int dimension9 = 90;
    cout << "Dimension Set 9: " << dimension9 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set9(dimension9, vector<double>(dimension9));
    vector<vector<double>> pod_matrix2_set9(dimension9, vector<double>(dimension9));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension9; i++)
    {
        for (int j = 0; j < dimension9; j++)
        {
            pod_matrix1_set9[i][j] = filler;
            pod_matrix2_set9[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set9);
    print_partial_matrix(pod_matrix2_set9);

    // --------------- FOURTH SET ------------------
    int dimension10 = 100;
    cout << "Dimension Set 10: " << dimension10 << endl
        << endl;

    vector<vector<double>> pod_matrix1_set10(dimension10, vector<double>(dimension10));
    vector<vector<double>> pod_matrix2_set10(dimension10, vector<double>(dimension10));

    // Fill input matrices
    filler = 1.0;
    // Set 1
    for (int i = 0; i < dimension10; i++)
    {
        for (int j = 0; j < dimension10; j++)
        {
            pod_matrix1_set10[i][j] = filler;
            pod_matrix2_set10[i][j] = static_cast<double>((j % 2) + 1);
            filler++;
        }
    }
    print_partial_matrix(pod_matrix1_set10);
    print_partial_matrix(pod_matrix2_set10);


    // Get all diagonals
    // Set 1
    vector<vector<double>> all_diagonal1_set1(dimension1, vector<double>(dimension1));
    vector<vector<double>> all_diagonal2_set1(dimension1, vector<double>(dimension1));

    for (int i = 0; i < dimension1; i++)
    {
        all_diagonal1_set1[i] = get_diagonal(i, pod_matrix1_set1);
        all_diagonal2_set1[i] = get_diagonal(i, pod_matrix2_set1);
    }

    cout << "Diagonal 1 Set 1 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set1);

    cout << "Diagonal 2 Set 1 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set1);

    // Set 2
    vector<vector<double>> all_diagonal1_set2(dimension2, vector<double>(dimension2));
    vector<vector<double>> all_diagonal2_set2(dimension2, vector<double>(dimension2));

    for (int i = 0; i < dimension2; i++)
    {
        all_diagonal1_set2[i] = get_diagonal(i, pod_matrix1_set2);
        all_diagonal2_set2[i] = get_diagonal(i, pod_matrix2_set2);
    }

    cout << "Diagonal 1 Set 2 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set1);

    cout << "Diagonal 2 Set 2 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set2);

    // Set 3
    vector<vector<double>> all_diagonal1_set3(dimension3, vector<double>(dimension3));
    vector<vector<double>> all_diagonal2_set3(dimension3, vector<double>(dimension3));

    for (int i = 0; i < dimension3; i++)
    {
        all_diagonal1_set3[i] = get_diagonal(i, pod_matrix1_set3);
        all_diagonal2_set3[i] = get_diagonal(i, pod_matrix2_set3);
    }

    cout << "Diagonal 1 Set 3 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set3);

    cout << "Diagonal 2 Set 3 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set3);

    // Set 4
    vector<vector<double>> all_diagonal1_set4(dimension4, vector<double>(dimension4));
    vector<vector<double>> all_diagonal2_set4(dimension4, vector<double>(dimension4));

    for (int i = 0; i < dimension4; i++)
    {
        all_diagonal1_set4[i] = get_diagonal(i, pod_matrix1_set4);
        all_diagonal2_set4[i] = get_diagonal(i, pod_matrix2_set4);
    }

    cout << "Diagonal 1 Set 4 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set4);

    cout << "Diagonal 2 Set 4 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set4);

    // Set 5
    vector<vector<double>> all_diagonal1_set5(dimension5, vector<double>(dimension5));
    vector<vector<double>> all_diagonal2_set5(dimension5, vector<double>(dimension5));

    for (int i = 0; i < dimension5; i++)
    {
        all_diagonal1_set5[i] = get_diagonal(i, pod_matrix1_set5);
        all_diagonal2_set5[i] = get_diagonal(i, pod_matrix2_set5);
    }

    cout << "Diagonal 1 Set 5 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set5);

    cout << "Diagonal 2 Set 5 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set5);

    // Set 6
    vector<vector<double>> all_diagonal1_set6(dimension6, vector<double>(dimension6));
    vector<vector<double>> all_diagonal2_set6(dimension6, vector<double>(dimension6));

    for (int i = 0; i < dimension6; i++)
    {
        all_diagonal1_set6[i] = get_diagonal(i, pod_matrix1_set6);
        all_diagonal2_set6[i] = get_diagonal(i, pod_matrix2_set6);
    }

    cout << "Diagonal 1 Set 6 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set6);

    cout << "Diagonal 2 Set 6 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set6);

    // Set 7
    vector<vector<double>> all_diagonal1_set7(dimension7, vector<double>(dimension7));
    vector<vector<double>> all_diagonal2_set7(dimension7, vector<double>(dimension7));

    for (int i = 0; i < dimension7; i++)
    {
        all_diagonal1_set7[i] = get_diagonal(i, pod_matrix1_set7);
        all_diagonal2_set7[i] = get_diagonal(i, pod_matrix2_set7);
    }

    cout << "Diagonal 1 Set 7 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set7);

    cout << "Diagonal 2 Set 7 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set7);

    // Set 8
    vector<vector<double>> all_diagonal1_set8(dimension8, vector<double>(dimension8));
    vector<vector<double>> all_diagonal2_set8(dimension8, vector<double>(dimension8));

    for (int i = 0; i < dimension8; i++)
    {
        all_diagonal1_set8[i] = get_diagonal(i, pod_matrix1_set8);
        all_diagonal2_set8[i] = get_diagonal(i, pod_matrix2_set8);
    }

    cout << "Diagonal 1 Set 8 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set8);

    cout << "Diagonal 2 Set 8 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set8);

    // Set 9
    vector<vector<double>> all_diagonal1_set9(dimension9, vector<double>(dimension9));
    vector<vector<double>> all_diagonal2_set9(dimension9, vector<double>(dimension9));

    for (int i = 0; i < dimension9; i++)
    {
        all_diagonal1_set9[i] = get_diagonal(i, pod_matrix1_set9);
        all_diagonal2_set9[i] = get_diagonal(i, pod_matrix2_set9);
    }

    cout << "Diagonal 1 Set 9 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set9);

    cout << "Diagonal 2 Set 9 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set9);

    // Set 10
    vector<vector<double>> all_diagonal1_set10(dimension10, vector<double>(dimension10));
    vector<vector<double>> all_diagonal2_set10(dimension10, vector<double>(dimension10));

    for (int i = 0; i < dimension10; i++)
    {
        all_diagonal1_set10[i] = get_diagonal(i, pod_matrix1_set10);
        all_diagonal2_set10[i] = get_diagonal(i, pod_matrix2_set10);
    }

    cout << "Diagonal 1 Set 10 Expected:" << endl;
    print_partial_matrix(all_diagonal1_set10);

    cout << "Diagonal 2 Set 10 Expected:" << endl;
    print_partial_matrix(all_diagonal2_set10);

    // Encode Matrices into vectors with Diagonals
    // Set 1
    vector<Plaintext> plain_matrix1_set1(dimension1), plain_matrix2_set1(dimension1);
    vector<Plaintext> plain_diagonal1_set1(dimension1), plain_diagonal2_set1(dimension1);

    for (int i = 0; i < dimension1; i++)
    {
        ckks_encoder.encode(pod_matrix1_set1[i], scale, plain_matrix1_set1[i]);
        ckks_encoder.encode(pod_matrix2_set1[i], scale, plain_matrix2_set1[i]);
        ckks_encoder.encode(all_diagonal1_set1[i], scale, plain_diagonal1_set1[i]);
        ckks_encoder.encode(all_diagonal2_set1[i], scale, plain_diagonal2_set1[i]);
    }

    cout << "Encoding Set 1 is Complete" << endl;

    // Set 2
    vector<Plaintext> plain_matrix1_set2(dimension2), plain_matrix2_set2(dimension2);
    vector<Plaintext> plain_diagonal1_set2(dimension2), plain_diagonal2_set2(dimension2);

    for (int i = 0; i < dimension2; i++)
    {
        ckks_encoder.encode(pod_matrix1_set2[i], scale, plain_matrix1_set2[i]);
        ckks_encoder.encode(pod_matrix2_set2[i], scale, plain_matrix2_set2[i]);
        ckks_encoder.encode(all_diagonal1_set2[i], scale, plain_diagonal1_set2[i]);
        ckks_encoder.encode(all_diagonal2_set2[i], scale, plain_diagonal2_set2[i]);
    }

    cout << "Encoding Set 2 is Complete" << endl;

    // Set 3
    vector<Plaintext> plain_matrix1_set3(dimension3), plain_matrix2_set3(dimension3);
    vector<Plaintext> plain_diagonal1_set3(dimension3), plain_diagonal2_set3(dimension3);

    for (int i = 0; i < dimension3; i++)
    {
        ckks_encoder.encode(pod_matrix1_set3[i], scale, plain_matrix1_set3[i]);
        ckks_encoder.encode(pod_matrix2_set3[i], scale, plain_matrix2_set3[i]);
        ckks_encoder.encode(all_diagonal1_set3[i], scale, plain_diagonal1_set3[i]);
        ckks_encoder.encode(all_diagonal2_set3[i], scale, plain_diagonal2_set3[i]);
    }

    cout << "Encoding Set 3 is Complete" << endl;

    // Set 4
    vector<Plaintext> plain_matrix1_set4(dimension4), plain_matrix2_set4(dimension4);
    vector<Plaintext> plain_diagonal1_set4(dimension4), plain_diagonal2_set4(dimension4);

    for (int i = 0; i < dimension4; i++)
    {
        ckks_encoder.encode(pod_matrix1_set4[i], scale, plain_matrix1_set4[i]);
        ckks_encoder.encode(pod_matrix2_set4[i], scale, plain_matrix2_set4[i]);
        ckks_encoder.encode(all_diagonal1_set4[i], scale, plain_diagonal1_set4[i]);
        ckks_encoder.encode(all_diagonal2_set4[i], scale, plain_diagonal2_set4[i]);
    }

    cout << "Encoding Set 4 is Complete" << endl;

    // Set 5
    vector<Plaintext> plain_matrix1_set5(dimension5), plain_matrix2_set5(dimension5);
    vector<Plaintext> plain_diagonal1_set5(dimension5), plain_diagonal2_set5(dimension5);

    for (int i = 0; i < dimension5; i++)
    {
        ckks_encoder.encode(pod_matrix1_set5[i], scale, plain_matrix1_set5[i]);
        ckks_encoder.encode(pod_matrix2_set5[i], scale, plain_matrix2_set5[i]);
        ckks_encoder.encode(all_diagonal1_set5[i], scale, plain_diagonal1_set5[i]);
        ckks_encoder.encode(all_diagonal2_set5[i], scale, plain_diagonal2_set5[i]);
    }

    cout << "Encoding Set 5 is Complete" << endl;

    // Set 6
    vector<Plaintext> plain_matrix1_set6(dimension6), plain_matrix2_set6(dimension6);
    vector<Plaintext> plain_diagonal1_set6(dimension6), plain_diagonal2_set6(dimension6);

    for (int i = 0; i < dimension6; i++)
    {
        ckks_encoder.encode(pod_matrix1_set6[i], scale, plain_matrix1_set6[i]);
        ckks_encoder.encode(pod_matrix2_set6[i], scale, plain_matrix2_set6[i]);
        ckks_encoder.encode(all_diagonal1_set6[i], scale, plain_diagonal1_set6[i]);
        ckks_encoder.encode(all_diagonal2_set6[i], scale, plain_diagonal2_set6[i]);
    }

    cout << "Encoding Set 6 is Complete" << endl;
    // Set 7
    vector<Plaintext> plain_matrix1_set7(dimension7), plain_matrix2_set7(dimension7);
    vector<Plaintext> plain_diagonal1_set7(dimension7), plain_diagonal2_set7(dimension7);

    for (int i = 0; i < dimension7; i++)
    {
        ckks_encoder.encode(pod_matrix1_set7[i], scale, plain_matrix1_set7[i]);
        ckks_encoder.encode(pod_matrix2_set7[i], scale, plain_matrix2_set7[i]);
        ckks_encoder.encode(all_diagonal1_set7[i], scale, plain_diagonal1_set7[i]);
        ckks_encoder.encode(all_diagonal2_set7[i], scale, plain_diagonal2_set7[i]);
    }

    cout << "Encoding Set 7 is Complete" << endl;

    // Set 8
    vector<Plaintext> plain_matrix1_set8(dimension8), plain_matrix2_set8(dimension8);
    vector<Plaintext> plain_diagonal1_set8(dimension8), plain_diagonal2_set8(dimension8);

    for (int i = 0; i < dimension8; i++)
    {
        ckks_encoder.encode(pod_matrix1_set8[i], scale, plain_matrix1_set8[i]);
        ckks_encoder.encode(pod_matrix2_set8[i], scale, plain_matrix2_set8[i]);
        ckks_encoder.encode(all_diagonal1_set8[i], scale, plain_diagonal1_set8[i]);
        ckks_encoder.encode(all_diagonal2_set8[i], scale, plain_diagonal2_set8[i]);
    }

    cout << "Encoding Set 8 is Complete" << endl;

    // Set 9
    vector<Plaintext> plain_matrix1_set9(dimension9), plain_matrix2_set9(dimension9);
    vector<Plaintext> plain_diagonal1_set9(dimension9), plain_diagonal2_set9(dimension9);

    for (int i = 0; i < dimension9; i++)
    {
        ckks_encoder.encode(pod_matrix1_set9[i], scale, plain_matrix1_set9[i]);
        ckks_encoder.encode(pod_matrix2_set9[i], scale, plain_matrix2_set9[i]);
        ckks_encoder.encode(all_diagonal1_set9[i], scale, plain_diagonal1_set9[i]);
        ckks_encoder.encode(all_diagonal2_set9[i], scale, plain_diagonal2_set9[i]);
    }

    cout << "Encoding Set 9 is Complete" << endl;

    // Set 10
    vector<Plaintext> plain_matrix1_set10(dimension10), plain_matrix2_set10(dimension10);
    vector<Plaintext> plain_diagonal1_set10(dimension10), plain_diagonal2_set10(dimension10);

    for (int i = 0; i < dimension10; i++)
    {
        ckks_encoder.encode(pod_matrix1_set10[i], scale, plain_matrix1_set10[i]);
        ckks_encoder.encode(pod_matrix2_set10[i], scale, plain_matrix2_set10[i]);
        ckks_encoder.encode(all_diagonal1_set10[i], scale, plain_diagonal1_set10[i]);
        ckks_encoder.encode(all_diagonal2_set10[i], scale, plain_diagonal2_set10[i]);
    }

    cout << "Encoding Set 10 is Complete" << endl;


    // Encrypt the matrices with Diagonals
    // Set 1
    vector<Ciphertext> cipher_matrix1_set1(dimension1), cipher_matrix2_set1(dimension1);
    vector<Ciphertext> cipher_diagonal1_set1(dimension1), cipher_diagonal2_set1(dimension1);

    for (unsigned int i = 0; i < dimension1; i++)
    {
        encryptor.encrypt(plain_matrix1_set1[i], cipher_matrix1_set1[i]);
        encryptor.encrypt(plain_matrix2_set1[i], cipher_matrix2_set1[i]);
        encryptor.encrypt(plain_diagonal1_set1[i], cipher_diagonal1_set1[i]);
        encryptor.encrypt(plain_diagonal2_set1[i], cipher_diagonal2_set1[i]);
    }
    cout << "Encrypting Set 1 is Complete" << endl;

    // Set 2
    vector<Ciphertext> cipher_matrix1_set2(dimension2), cipher_matrix2_set2(dimension2);
    vector<Ciphertext> cipher_diagonal1_set2(dimension2), cipher_diagonal2_set2(dimension2);

    for (unsigned int i = 0; i < dimension2; i++)
    {
        encryptor.encrypt(plain_matrix1_set2[i], cipher_matrix1_set2[i]);
        encryptor.encrypt(plain_matrix2_set2[i], cipher_matrix2_set2[i]);
        encryptor.encrypt(plain_diagonal1_set2[i], cipher_diagonal1_set2[i]);
        encryptor.encrypt(plain_diagonal2_set2[i], cipher_diagonal2_set2[i]);
    }
    cout << "Encrypting Set 2 is Complete" << endl;

    // Set 3
    vector<Ciphertext> cipher_matrix1_set3(dimension3), cipher_matrix2_set3(dimension3);
    vector<Ciphertext> cipher_diagonal1_set3(dimension3), cipher_diagonal2_set3(dimension3);

    for (unsigned int i = 0; i < dimension3; i++)
    {
        encryptor.encrypt(plain_matrix1_set3[i], cipher_matrix1_set3[i]);
        encryptor.encrypt(plain_matrix2_set3[i], cipher_matrix2_set3[i]);
        encryptor.encrypt(plain_diagonal1_set3[i], cipher_diagonal1_set3[i]);
        encryptor.encrypt(plain_diagonal2_set3[i], cipher_diagonal2_set3[i]);
    }
    cout << "Encrypting Set 3 is Complete" << endl;

    // Set 4
    vector<Ciphertext> cipher_matrix1_set4(dimension4), cipher_matrix2_set4(dimension4);
    vector<Ciphertext> cipher_diagonal1_set4(dimension4), cipher_diagonal2_set4(dimension4);

    for (unsigned int i = 0; i < dimension4; i++)
    {
        encryptor.encrypt(plain_matrix1_set4[i], cipher_matrix1_set4[i]);
        encryptor.encrypt(plain_matrix2_set4[i], cipher_matrix2_set4[i]);
        encryptor.encrypt(plain_diagonal1_set4[i], cipher_diagonal1_set4[i]);
        encryptor.encrypt(plain_diagonal2_set4[i], cipher_diagonal2_set4[i]);
    }
    cout << "Encrypting Set 4 is Complete" << endl;

    // Set 5
    vector<Ciphertext> cipher_matrix1_set5(dimension5), cipher_matrix2_set5(dimension5);
    vector<Ciphertext> cipher_diagonal1_set5(dimension5), cipher_diagonal2_set5(dimension5);

    for (unsigned int i = 0; i < dimension5; i++)
    {
        encryptor.encrypt(plain_matrix1_set5[i], cipher_matrix1_set5[i]);
        encryptor.encrypt(plain_matrix2_set5[i], cipher_matrix2_set5[i]);
        encryptor.encrypt(plain_diagonal1_set5[i], cipher_diagonal1_set5[i]);
        encryptor.encrypt(plain_diagonal2_set5[i], cipher_diagonal2_set5[i]);
    }
    cout << "Encrypting Set 5 is Complete" << endl;

    // Set 6
    vector<Ciphertext> cipher_matrix1_set6(dimension6), cipher_matrix2_set6(dimension6);
    vector<Ciphertext> cipher_diagonal1_set6(dimension6), cipher_diagonal2_set6(dimension6);

    for (unsigned int i = 0; i < dimension6; i++)
    {
        encryptor.encrypt(plain_matrix1_set6[i], cipher_matrix1_set6[i]);
        encryptor.encrypt(plain_matrix2_set6[i], cipher_matrix2_set6[i]);
        encryptor.encrypt(plain_diagonal1_set6[i], cipher_diagonal1_set6[i]);
        encryptor.encrypt(plain_diagonal2_set6[i], cipher_diagonal2_set6[i]);
    }
    cout << "Encrypting Set 6 is Complete" << endl;

    // Set 7
    vector<Ciphertext> cipher_matrix1_set7(dimension7), cipher_matrix2_set7(dimension7);
    vector<Ciphertext> cipher_diagonal1_set7(dimension7), cipher_diagonal2_set7(dimension7);

    for (unsigned int i = 0; i < dimension7; i++)
    {
        encryptor.encrypt(plain_matrix1_set7[i], cipher_matrix1_set7[i]);
        encryptor.encrypt(plain_matrix2_set7[i], cipher_matrix2_set7[i]);
        encryptor.encrypt(plain_diagonal1_set7[i], cipher_diagonal1_set7[i]);
        encryptor.encrypt(plain_diagonal2_set7[i], cipher_diagonal2_set7[i]);
    }
    cout << "Encrypting Set 7 is Complete" << endl;

    // Set 8
    vector<Ciphertext> cipher_matrix1_set8(dimension8), cipher_matrix2_set8(dimension8);
    vector<Ciphertext> cipher_diagonal1_set8(dimension8), cipher_diagonal2_set8(dimension8);

    for (unsigned int i = 0; i < dimension8; i++)
    {
        encryptor.encrypt(plain_matrix1_set8[i], cipher_matrix1_set8[i]);
        encryptor.encrypt(plain_matrix2_set8[i], cipher_matrix2_set8[i]);
        encryptor.encrypt(plain_diagonal1_set8[i], cipher_diagonal1_set8[i]);
        encryptor.encrypt(plain_diagonal2_set8[i], cipher_diagonal2_set8[i]);
    }
    cout << "Encrypting Set 8 is Complete" << endl;

    // Set 9
    vector<Ciphertext> cipher_matrix1_set9(dimension9), cipher_matrix2_set9(dimension9);
    vector<Ciphertext> cipher_diagonal1_set9(dimension9), cipher_diagonal2_set9(dimension9);

    for (unsigned int i = 0; i < dimension9; i++)
    {
        encryptor.encrypt(plain_matrix1_set9[i], cipher_matrix1_set9[i]);
        encryptor.encrypt(plain_matrix2_set9[i], cipher_matrix2_set9[i]);
        encryptor.encrypt(plain_diagonal1_set9[i], cipher_diagonal1_set9[i]);
        encryptor.encrypt(plain_diagonal2_set9[i], cipher_diagonal2_set9[i]);
    }
    cout << "Encrypting Set 9 is Complete" << endl;

    // Set 10
    vector<Ciphertext> cipher_matrix1_set10(dimension10), cipher_matrix2_set10(dimension10);
    vector<Ciphertext> cipher_diagonal1_set10(dimension10), cipher_diagonal2_set10(dimension10);

    for (unsigned int i = 0; i < dimension10; i++)
    {
        encryptor.encrypt(plain_matrix1_set10[i], cipher_matrix1_set10[i]);
        encryptor.encrypt(plain_matrix2_set10[i], cipher_matrix2_set10[i]);
        encryptor.encrypt(plain_diagonal1_set10[i], cipher_diagonal1_set10[i]);
        encryptor.encrypt(plain_diagonal2_set10[i], cipher_diagonal2_set10[i]);
    }
    cout << "Encrypting Set 10 is Complete" << endl;

    // ------------- FIRST COMPUTATION ----------------
    outf << "# index 0" << endl;
    outf << "# C_Vec . P_Mat" << endl;

    /*
        // test decrypt here
        for (unsigned int i = 0; i < dimension1; i++)
        {
            decryptor.decrypt(cipher_matrix1_set1[i], plain_matrix1_set1[i]);
            decryptor.decrypt(cipher_matrix2_set1[i], plain_matrix2_set1[i]);
            decryptor.decrypt(cipher_diagonal1_set1[i], plain_diagonal1_set1[i]);
            decryptor.decrypt(cipher_diagonal2_set1[i], plain_diagonal2_set1[i]);
        }
        // test decode here
        // test decrypt here
        for (unsigned int i = 0; i < dimension1; i++)
        {
            ckks_encoder.decode(plain_diagonal1_set1[i], all_diagonal_1[i]);
            ckks_encoder.decode(plain_diagonal2_set1[i], all_diagonal_2[i]);
        }
        // test print output
        cout << "\nDiagonal Set 1 Result:" << endl;
        for (unsigned int i = 0; i < dimension1; i++)
        {
            cout << "\t[";
            for (unsigned int j = 0; j < dimension1 - 1; j++)
            {
                cout << all_diagonal_1[i][j] << ", ";
            }
            cout << all_diagonal_1[i][dimension1 - 1];
            cout << "]" << endl;
        }
        cout << "\n"
            << endl;
    */
    // Test LinearTransform here
    // Set 1
    auto start_comp1_set1 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime1_set1 = Linear_Transform_Plain(cipher_matrix1_set1[0], plain_diagonal1_set1, gal_keys, params);
    auto stop_comp1_set1 = chrono::high_resolution_clock::now();

    /*
        UNCOMMENT TO DEBUG LINEAR TRANSFORM FUNCTION
        // Fill ct
        Ciphertext ct_rotated;
        evaluator.rotate_vector(cipher_matrix1_set1[0], -dimension1, gal_keys, ct_rotated);
        Ciphertext ct;
        evaluator.add(cipher_matrix1_set1[0], ct_rotated, ct);
        // Add epsilon to avoid negative numbers
        vector<double> epsilon_vec(poly_modulus_degree / 2);
        for (int i = 0; i < epsilon_vec.size(); i++)
        {
            epsilon_vec[i] = 0.0000;
        }
        Plaintext epsilon_plain;
        ckks_encoder.encode(epsilon_vec, scale, epsilon_plain);
        evaluator.add_plain_inplace(ct, epsilon_plain);
        // test fill ct
        Plaintext test_fill;
        decryptor.decrypt(ct, test_fill);
        vector<double> out_fill;
        ckks_encoder.decode(test_fill, out_fill);
        cout << "Filled CT:\n"
            << endl;
        for (int i = 0; i < dimension1; i++)
        {
            cout << out_fill[i] << ", ";
        }
        cout << "\n"
            << endl;
        Ciphertext ct_prime;
        // ct` = CMult(ct, u0)
        evaluator.multiply_plain(ct, plain_diagonal1_set1[0], ct_prime);
        // test mult plain 0
        Plaintext test_0;
        decryptor.decrypt(ct_prime, test_0);
        vector<double> out_test_0;
        ckks_encoder.decode(test_0, out_test_0);
        cout << "CT_Prime 0 :\n"
            << endl;
        for (int i = 0; i < dimension1; i++)
        {
            cout << out_test_0[i] << ", ";
        }
        cout << "\n"
            << endl;
        for (int l = 1; l < dimension1; l++)
        {
            // ct` = Add(ct`, CMult(Rot(ct, l), ul))
            Ciphertext temp_rot;
            Ciphertext temp_mul;
            evaluator.rotate_vector(ct, l, gal_keys, temp_rot);
            evaluator.multiply_plain(temp_rot, plain_diagonal1_set1[l], temp_mul);
            evaluator.add_inplace(ct_prime, temp_mul);
            // test decrypt
            Plaintext temp_rot_plain;
            Plaintext temp_mul_plain;
            Plaintext temp_ct_prime;
            decryptor.decrypt(temp_rot, temp_rot_plain);
            decryptor.decrypt(temp_mul, temp_mul_plain);
            decryptor.decrypt(ct_prime, temp_ct_prime);
            // test decode
            vector<double> test_out_rot, test_out_mul, test_ct_prime;
            vector<double> test_diag;
            ckks_encoder.decode(temp_ct_prime, test_ct_prime);
            ckks_encoder.decode(temp_mul_plain, test_out_mul);
            ckks_encoder.decode(temp_rot_plain, test_out_rot);
            ckks_encoder.decode(plain_diagonal1_set1[l], test_diag);
            cout << "Rotation " << l << "\n"
                << endl;
            cout << "\nrotated vec:\n\t[";
            for (int j = 0; j < dimension1; j++)
            {
                cout << test_out_rot[j] << ", ";
            }
            cout << "\nDiagonal vec:\n\t[";
            for (int j = 0; j < dimension1; j++)
            {
                cout << test_diag[j] << ", ";
            }
            cout << "\nMult vec vec:\n\t[";
            for (int j = 0; j < dimension1; j++)
            {
                cout << test_out_mul[j] << ", ";
            }
            cout << "\nCt_prime vec:\n\t[";
            for (int j = 0; j < dimension1; j++)
            {
                cout << test_ct_prime[j] << ", ";
            }
            cout << "\n"
                << endl;
        }
    */
    // Decrypt
    Plaintext pt_result1_set1;
    decryptor.decrypt(ct_prime1_set1, pt_result1_set1);

    // Decode
    vector<double> output_result1_set1;
    ckks_encoder.decode(pt_result1_set1, output_result1_set1);

    auto duration_comp1_set1 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set1 - start_comp1_set1);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set1.count() << " microseconds" << endl;
    outf << "10\t\t" << duration_comp1_set1.count() << endl;

    cout << "Linear Transformation Set 1 Result:" << endl;
    print_partial_vector(output_result1_set1, dimension1);

    // Check result
    cout << "Expected output Set 1: " << endl;

    test_Linear_Transformation(dimension1, pod_matrix1_set1, pod_matrix1_set1[0]);

    // Set 2
    auto start_comp1_set2 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set2 = Linear_Transform_Plain(cipher_matrix1_set2[0], plain_diagonal1_set2, gal_keys, params);
    auto stop_comp1_set2 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result_set2;
    decryptor.decrypt(ct_prime_set2, pt_result_set2);

    // Decode
    vector<double> output_result_set2;
    ckks_encoder.decode(pt_result_set2, output_result_set2);

    auto duration_comp1_set2 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set2 - start_comp1_set2);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set2.count() << " microseconds" << endl;
    outf << "100\t\t" << duration_comp1_set2.count() << endl;

    cout << "Linear Transformation Set 2 Result:" << endl;
    print_partial_vector(output_result_set2, dimension2);

    // Check result
    cout << "Expected output Set 2: " << endl;

    test_Linear_Transformation(dimension2, pod_matrix1_set2, pod_matrix1_set2[0]);

    // Set 3
    auto start_comp1_set3 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set3 = Linear_Transform_Plain(cipher_matrix1_set3[0], plain_diagonal1_set3, gal_keys, params);
    auto stop_comp1_set3 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set3;
    decryptor.decrypt(ct_prime_set3, pt_result_set3);

    // Decode
    vector<double> output_result_set3;
    ckks_encoder.decode(pt_result_set3, output_result_set3);

    auto duration_comp1_set3 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set3 - start_comp1_set3);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set3.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set3.count() << endl;

    cout << "Linear Transformation Set 3 Result:" << endl;
    print_partial_vector(output_result_set3, dimension3);

    // Check result
    cout << "Expected output Set 3: " << endl;
    test_Linear_Transformation(dimension3, pod_matrix1_set3, pod_matrix1_set3[0]);

    // Set 4
    auto start_comp1_set4 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set4 = Linear_Transform_Plain(cipher_matrix1_set4[0], plain_diagonal1_set4, gal_keys, params);
    auto stop_comp1_set4 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set4;
    decryptor.decrypt(ct_prime_set4, pt_result_set4);

    // Decode
    vector<double> output_result_set4;
    ckks_encoder.decode(pt_result_set4, output_result_set4);

    auto duration_comp1_set4 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set4 - start_comp1_set4);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set4.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set4.count() << endl;

    cout << "Linear Transformation Set 4 Result:" << endl;
    print_partial_vector(output_result_set4, dimension4);

    // Check result
    cout << "Expected output Set 4: " << endl;
    test_Linear_Transformation(dimension4, pod_matrix1_set4, pod_matrix1_set4[0]);

    // Set 5
    auto start_comp1_set5 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set5 = Linear_Transform_Plain(cipher_matrix1_set5[0], plain_diagonal1_set5, gal_keys, params);
    auto stop_comp1_set5 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set5;
    decryptor.decrypt(ct_prime_set5, pt_result_set5);

    // Decode
    vector<double> output_result_set5;
    ckks_encoder.decode(pt_result_set5, output_result_set5);

    auto duration_comp1_set5 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set5 - start_comp1_set5);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set5.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set5.count() << endl;

    cout << "Linear Transformation Set 5 Result:" << endl;
    print_partial_vector(output_result_set5, dimension5);

    // Check result
    cout << "Expected output Set 5: " << endl;
    test_Linear_Transformation(dimension5, pod_matrix1_set5, pod_matrix1_set5[0]);

    // Set 6
    auto start_comp1_set6 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set6 = Linear_Transform_Plain(cipher_matrix1_set6[0], plain_diagonal1_set6, gal_keys, params);
    auto stop_comp1_set6 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set6;
    decryptor.decrypt(ct_prime_set6, pt_result_set6);

    // Decode
    vector<double> output_result_set6;
    ckks_encoder.decode(pt_result_set6, output_result_set6);

    auto duration_comp1_set6 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set6 - start_comp1_set6);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set6.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set6.count() << endl;

    cout << "Linear Transformation Set 6 Result:" << endl;
    print_partial_vector(output_result_set6, dimension6);

    // Check result
    cout << "Expected output Set 6: " << endl;
    test_Linear_Transformation(dimension6, pod_matrix1_set6, pod_matrix1_set6[0]);

    // Set 7
    auto start_comp1_set7 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set7 = Linear_Transform_Plain(cipher_matrix1_set7[0], plain_diagonal1_set7, gal_keys, params);
    auto stop_comp1_set7 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set7;
    decryptor.decrypt(ct_prime_set7, pt_result_set7);

    // Decode
    vector<double> output_result_set7;
    ckks_encoder.decode(pt_result_set7, output_result_set7);

    auto duration_comp1_set7 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set7 - start_comp1_set7);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set7.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set7.count() << endl;

    cout << "Linear Transformation Set 7 Result:" << endl;
    print_partial_vector(output_result_set7, dimension7);

    // Check result
    cout << "Expected output Set 7: " << endl;
    test_Linear_Transformation(dimension7, pod_matrix1_set7, pod_matrix1_set7[0]);

    // Set 8
    auto start_comp1_set8 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set8 = Linear_Transform_Plain(cipher_matrix1_set8[0], plain_diagonal1_set8, gal_keys, params);
    auto stop_comp1_set8 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set8;
    decryptor.decrypt(ct_prime_set8, pt_result_set8);

    // Decode
    vector<double> output_result_set8;
    ckks_encoder.decode(pt_result_set8, output_result_set8);

    auto duration_comp1_set8 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set8 - start_comp1_set8);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set8.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set8.count() << endl;

    cout << "Linear Transformation Set 8 Result:" << endl;
    print_partial_vector(output_result_set8, dimension8);

    // Check result
    cout << "Expected output Set 8: " << endl;
    test_Linear_Transformation(dimension8, pod_matrix1_set8, pod_matrix1_set8[0]);

    // Set 9
    auto start_comp1_set9 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set9 = Linear_Transform_Plain(cipher_matrix1_set9[0], plain_diagonal1_set9, gal_keys, params);
    auto stop_comp1_set9 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set9;
    decryptor.decrypt(ct_prime_set9, pt_result_set9);

    // Decode
    vector<double> output_result_set9;
    ckks_encoder.decode(pt_result_set9, output_result_set9);

    auto duration_comp1_set9 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set9 - start_comp1_set9);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set9.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set9.count() << endl;

    cout << "Linear Transformation Set 9 Result:" << endl;
    print_partial_vector(output_result_set9, dimension9);

    // Check result
    cout << "Expected output Set 9: " << endl;
    test_Linear_Transformation(dimension9, pod_matrix1_set9, pod_matrix1_set9[0]);

    // Set 10
    auto start_comp1_set10 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime_set10 = Linear_Transform_Plain(cipher_matrix1_set10[0], plain_diagonal1_set10, gal_keys, params);
    auto stop_comp1_set10 = chrono::high_resolution_clock::now();
    // Decrypt
    Plaintext pt_result_set10;
    decryptor.decrypt(ct_prime_set10, pt_result_set10);

    // Decode
    vector<double> output_result_set10;
    ckks_encoder.decode(pt_result_set10, output_result_set10);

    auto duration_comp1_set10 = chrono::duration_cast<chrono::microseconds>(stop_comp1_set10 - start_comp1_set10);
    cout << "\nTime to compute C_vec . P_mat: " << duration_comp1_set10.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp1_set10.count() << endl;

    cout << "Linear Transformation Set 10 Result:" << endl;
    print_partial_vector(output_result_set10, dimension10);

    // Check result
    cout << "Expected output Set 10: " << endl;
    test_Linear_Transformation(dimension10, pod_matrix1_set10, pod_matrix1_set10[0]);

    outf << "\n"
        << endl;

    // ------------- TIME DIFFERENCE----------------

    cout << "\nTime to compute C_vec . P_mat: " << endl;

    cout << "\n10 - 20 : " << duration_comp1_set2.count() - duration_comp1_set1.count() << " microseconds" << endl;
    cout << "\n20 - 30 : " << duration_comp1_set3.count() - duration_comp1_set2.count() << " microseconds" << endl;
    cout << "\n30 - 40 : " << duration_comp1_set4.count() - duration_comp1_set3.count() << " microseconds" << endl;
    cout << "\n40 - 50 : " << duration_comp1_set5.count() - duration_comp1_set4.count() << " microseconds" << endl;
    cout << "\n50 - 60 : " << duration_comp1_set6.count() - duration_comp1_set5.count() << " microseconds" << endl;
    cout << "\n60 - 70 : " << duration_comp1_set7.count() - duration_comp1_set6.count() << " microseconds" << endl;
    cout << "\n70 - 80 : " << duration_comp1_set8.count() - duration_comp1_set7.count() << " microseconds" << endl;
    cout << "\n80 - 90 : " << duration_comp1_set9.count() - duration_comp1_set8.count() << " microseconds" << endl;
    cout << "\n90 - 100 : " << duration_comp1_set10.count() - duration_comp1_set9.count() << " microseconds" << endl;

    // ------------- SECOND COMPUTATION ----------------
    outf << "# index 1" << endl;
    outf << "# C_Vec . C_Mat" << endl;
    // Set 1
    auto start_comp2_set1 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set1 = Linear_Transform_Cipher(cipher_matrix1_set1[0], cipher_diagonal1_set1, gal_keys, params);
    auto stop_comp2_set1 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set1;
    decryptor.decrypt(ct_prime2_set1, pt_result2_set1);

    // Decode
    vector<double> output_result2_set1;
    ckks_encoder.decode(pt_result2_set1, output_result2_set1);

    auto duration_comp2_set1 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set1 - start_comp2_set1);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set1.count() << " microseconds" << endl;
    outf << "10\t\t" << duration_comp2_set1.count() << endl;

    cout << "Linear Transformation Set 1 Result:" << endl;
    print_partial_vector(output_result2_set1, dimension1);

    // Check result
    cout << "Expected output Set 1: " << endl;

    test_Linear_Transformation(dimension1, pod_matrix1_set1, pod_matrix1_set1[0]);

    // Set 2
    auto start_comp2_set2 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set2 = Linear_Transform_Cipher(cipher_matrix1_set2[0], cipher_diagonal1_set2, gal_keys, params);
    auto stop_comp2_set2 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set2;
    decryptor.decrypt(ct_prime2_set2, pt_result2_set2);

    // Decode
    vector<double> output_result2_set2;
    ckks_encoder.decode(pt_result2_set2, output_result2_set2);

    auto duration_comp2_set2 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set2 - start_comp2_set2);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set2.count() << " microseconds" << endl;
    outf << "100\t\t" << duration_comp2_set2.count() << endl;

    cout << "Linear Transformation Set 2 Result:" << endl;
    print_partial_vector(output_result2_set2, dimension2);

    // Check result
    cout << "Expected output Set 2: " << endl;

    test_Linear_Transformation(dimension2, pod_matrix1_set2, pod_matrix1_set2[0]);

    // Set 3
    auto start_comp2_set3 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set3 = Linear_Transform_Cipher(cipher_matrix1_set3[0], cipher_diagonal1_set3, gal_keys, params);
    auto stop_comp2_set3 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set3;
    decryptor.decrypt(ct_prime2_set3, pt_result2_set3);

    // Decode
    vector<double> output_result2_set3;
    ckks_encoder.decode(pt_result2_set3, output_result2_set3);

    auto duration_comp2_set3 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set3 - start_comp2_set3);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set3.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set3.count() << endl;

    cout << "Linear Transformation Set 3 Result:" << endl;
    print_partial_vector(output_result2_set3, dimension3);

    // Check result
    cout << "Expected output Set 3: " << endl;

    test_Linear_Transformation(dimension3, pod_matrix1_set3, pod_matrix1_set3[0]);

    // Set 4
    auto start_comp2_set4 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set4 = Linear_Transform_Cipher(cipher_matrix1_set4[0], cipher_diagonal1_set4, gal_keys, params);
    auto stop_comp2_set4 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set4;
    decryptor.decrypt(ct_prime2_set4, pt_result2_set4);

    // Decode
    vector<double> output_result2_set4;
    ckks_encoder.decode(pt_result2_set4, output_result2_set4);

    auto duration_comp2_set4 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set4 - start_comp2_set4);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set4.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set4.count() << endl;

    cout << "Linear Transformation Set 4 Result:" << endl;
    print_partial_vector(output_result2_set4, dimension4);

    // Check result
    cout << "Expected output Set 4: " << endl;

    test_Linear_Transformation(dimension4, pod_matrix1_set4, pod_matrix1_set4[0]);

    // Set 5
    auto start_comp2_set5 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set5 = Linear_Transform_Cipher(cipher_matrix1_set5[0], cipher_diagonal1_set5, gal_keys, params);
    auto stop_comp2_set5 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set5;
    decryptor.decrypt(ct_prime2_set5, pt_result2_set5);

    // Decode
    vector<double> output_result2_set5;
    ckks_encoder.decode(pt_result2_set5, output_result2_set5);

    auto duration_comp2_set5 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set5 - start_comp2_set5);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set5.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set5.count() << endl;

    cout << "Linear Transformation Set 5 Result:" << endl;
    print_partial_vector(output_result2_set5, dimension5);

    // Check result
    cout << "Expected output Set 5: " << endl;

    test_Linear_Transformation(dimension5, pod_matrix1_set5, pod_matrix1_set5[0]);

    // Set 6
    auto start_comp2_set6 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set6 = Linear_Transform_Cipher(cipher_matrix1_set6[0], cipher_diagonal1_set6, gal_keys, params);
    auto stop_comp2_set6 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set6;
    decryptor.decrypt(ct_prime2_set6, pt_result2_set6);

    // Decode
    vector<double> output_result2_set6;
    ckks_encoder.decode(pt_result2_set6, output_result2_set6);

    auto duration_comp2_set6 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set6 - start_comp2_set6);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set6.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set6.count() << endl;

    cout << "Linear Transformation Set 6 Result:" << endl;
    print_partial_vector(output_result2_set6, dimension6);

    // Check result
    cout << "Expected output Set 6: " << endl;

    test_Linear_Transformation(dimension6, pod_matrix1_set6, pod_matrix1_set6[0]);

    // Set 7
    auto start_comp2_set7 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set7 = Linear_Transform_Cipher(cipher_matrix1_set7[0], cipher_diagonal1_set7, gal_keys, params);
    auto stop_comp2_set7 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set7;
    decryptor.decrypt(ct_prime2_set7, pt_result2_set7);

    // Decode
    vector<double> output_result2_set7;
    ckks_encoder.decode(pt_result2_set7, output_result2_set7);

    auto duration_comp2_set7 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set7 - start_comp2_set7);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set7.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set7.count() << endl;

    cout << "Linear Transformation Set 7 Result:" << endl;
    print_partial_vector(output_result2_set7, dimension7);

    // Check result
    cout << "Expected output Set 7: " << endl;

    test_Linear_Transformation(dimension7, pod_matrix1_set7, pod_matrix1_set7[0]);

    // Set 8
    auto start_comp2_set8 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set8 = Linear_Transform_Cipher(cipher_matrix1_set8[0], cipher_diagonal1_set8, gal_keys, params);
    auto stop_comp2_set8 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set8;
    decryptor.decrypt(ct_prime2_set8, pt_result2_set8);

    // Decode
    vector<double> output_result2_set8;
    ckks_encoder.decode(pt_result2_set8, output_result2_set8);

    auto duration_comp2_set8 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set8 - start_comp2_set8);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set8.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set8.count() << endl;

    cout << "Linear Transformation Set 8 Result:" << endl;
    print_partial_vector(output_result2_set8, dimension8);

    // Check result
    cout << "Expected output Set 8: " << endl;

    test_Linear_Transformation(dimension8, pod_matrix1_set8, pod_matrix1_set8[0]);

    // Set 9
    auto start_comp2_set9 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set9 = Linear_Transform_Cipher(cipher_matrix1_set9[0], cipher_diagonal1_set9, gal_keys, params);
    auto stop_comp2_set9 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set9;
    decryptor.decrypt(ct_prime2_set9, pt_result2_set9);

    // Decode
    vector<double> output_result2_set9;
    ckks_encoder.decode(pt_result2_set9, output_result2_set9);

    auto duration_comp2_set9 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set9 - start_comp2_set9);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set9.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set9.count() << endl;

    cout << "Linear Transformation Set 9 Result:" << endl;
    print_partial_vector(output_result2_set9, dimension9);

    // Check result
    cout << "Expected output Set 9: " << endl;

    test_Linear_Transformation(dimension9, pod_matrix1_set9, pod_matrix1_set9[0]);

    // Set 10
    auto start_comp2_set10 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set10 = Linear_Transform_Cipher(cipher_matrix1_set10[0], cipher_diagonal1_set10, gal_keys, params);
    auto stop_comp2_set10 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext pt_result2_set10;
    decryptor.decrypt(ct_prime2_set10, pt_result2_set10);

    // Decode
    vector<double> output_result2_set10;
    ckks_encoder.decode(pt_result2_set10, output_result2_set10);

    auto duration_comp2_set10 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set10 - start_comp2_set10);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set10.count() << " microseconds" << endl;
    outf << "1000\t\t" << duration_comp2_set10.count() << endl;

    cout << "Linear Transformation Set 10 Result:" << endl;
    print_partial_vector(output_result2_set10, dimension10);

    // Check result
    cout << "Expected output Set 10: " << endl;

    test_Linear_Transformation(dimension10, pod_matrix1_set10, pod_matrix1_set10[0]);

    // ------------- TIME DIFFERENCE----------------

    cout << "\nTime to compute C_vec . C_mat: " << endl;

    cout << "\n10 - 20 : " << duration_comp2_set2.count() - duration_comp2_set1.count() << " microseconds" << endl;
    cout << "\n20 - 30 : " << duration_comp2_set3.count() - duration_comp2_set2.count() << " microseconds" << endl;
    cout << "\n30 - 40 : " << duration_comp2_set4.count() - duration_comp2_set3.count() << " microseconds" << endl;
    cout << "\n40 - 50 : " << duration_comp2_set5.count() - duration_comp2_set4.count() << " microseconds" << endl;
    cout << "\n50 - 60 : " << duration_comp2_set6.count() - duration_comp2_set5.count() << " microseconds" << endl;
    cout << "\n60 - 70 : " << duration_comp2_set7.count() - duration_comp2_set6.count() << " microseconds" << endl;
    cout << "\n70 - 80 : " << duration_comp2_set8.count() - duration_comp2_set7.count() << " microseconds" << endl;
    cout << "\n80 - 90 : " << duration_comp2_set9.count() - duration_comp2_set8.count() << " microseconds" << endl;
    cout << "\n90 - 100 : " << duration_comp2_set10.count() - duration_comp2_set9.count() << " microseconds" << endl;

    outf << "\n"
        << endl;
    outf.close();
}

int main()
{
    TimesForDifferentSets(8192);

    return 0;
}