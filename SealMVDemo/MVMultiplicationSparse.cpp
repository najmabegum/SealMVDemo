
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Common.h"
#include "seal/seal.h"


using namespace std;
using namespace seal;

void MVMultiplicationSparse(int dimension, int diagonalStart, int diagonalEnd, bool rescale)
{
    
    /*Set Seal context*/
    size_t poly_modulus_degree = 8192;
    EncryptionParameters params(scheme_type::ckks);
    params.set_poly_modulus_degree(poly_modulus_degree);
    cout << "MAX BIT COUNT: " << CoeffModulus::MaxBitCount(poly_modulus_degree) << endl;    
    params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    auto context = SEALContext::SEALContext(params);

    KeyGenerator keygen(context);
    auto sk = keygen.secret_key();
    PublicKey pk;
    keygen.create_public_key(pk);
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, pk);
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);

    // Create CKKS encoder
    CKKSEncoder ckks_encoder(context);
    // Create scale
    cout << "Coeff Modulus Back Value: " << params.coeff_modulus().back().value() << endl;
    double scale = pow(2.0, 40);


    // --------------- MATRIX SET ------------------
    
    cout << "Dimension of the Matrix " << dimension << endl
        << endl;

    vector<vector<double>> plain_matrix(dimension, vector<double>(dimension));

    // Fill input matrices
    double filler = 1.0;
    int counter = 0;

    // Set 1    
    for (int i = 0; i < dimension; i++)
    {
        counter = 0;
        for (int j = 0; j < dimension; j++)
        {
            if (counter < 20)
            {
                if (j >= diagonalStart && j <= diagonalEnd)
                {
                    plain_matrix[i][j] = filler;
                    filler++;
                    counter++;
                    if (j == diagonalEnd)
                    {
                        diagonalStart++;
                        diagonalEnd++;
                    }
                }
                else
                {
                    plain_matrix[i][j] = 0;
                }
            }            
            else
            {
                plain_matrix[i][j] = 0;
            }                        
        }
    }

    cout << "This is the input Matrix: " << endl;
    print_full_matrix(plain_matrix);
    
    // Get all diagonals
    vector<vector<double>> matrix_diagonals(dimension, vector<double>(dimension));

    for (int i = 0; i < dimension; i++)
    {
        matrix_diagonals[i] = get_diagonal(i, plain_matrix);
    }

    cout << "Matrix Diagonals: " << endl;
    print_partial_matrix(matrix_diagonals);


    // Encode Matrices into vectors with Diagonals 
    vector<Plaintext> encoded_diagonals(dimension);

    for (int i = 0; i < dimension; i++)
    {
        ckks_encoder.encode(matrix_diagonals[i], scale, encoded_diagonals[i]);
    }

    cout << "Encoding Complete" << endl;

    // Encrypt the matrices with Diagonals 
    vector<Ciphertext> encrypted_diagonals(dimension);

    for (unsigned int i = 0; i < dimension; i++)
    {

        encryptor.encrypt(encoded_diagonals[i], encrypted_diagonals[i]);

    }

    cout << "Encrypting Complete" << endl;

    cout << "Encoding and encrypting the Vector to be used for the multiplication" << endl;

    Plaintext encoded_vector;
    Ciphertext encrypted_vector;
    ckks_encoder.encode(plain_matrix[0], scale, encoded_vector);
    encryptor.encrypt(encoded_vector, encrypted_vector);

    cout << "Encoding and Encryption of the Vector complete.";

  

    // Linear Transform on Cipher
    auto start_comp2_set2 = chrono::high_resolution_clock::now();
    Ciphertext ct_prime2_set2 = Linear_Transform_Cipher(encrypted_vector, encrypted_diagonals, gal_keys, params, rescale, relin_keys);
    auto stop_comp2_set2 = chrono::high_resolution_clock::now();

    // Decrypt
    Plaintext encoded_result;
    decryptor.decrypt(ct_prime2_set2, encoded_result);

    // Decode
    vector<double> plain_result;
    ckks_encoder.decode(encoded_result, plain_result);

    auto duration_comp2_set2 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set2 - start_comp2_set2);
    cout << "\nTime to compute C_vec . C_mat: " << duration_comp2_set2.count() << " microseconds" << endl;

    cout << "HE Result:" << endl;
    print_partial_vector(plain_result, dimension);

    // Check result
    cout << "Expected Result: " << endl;

    test_Linear_Transformation(dimension, plain_matrix, plain_matrix[0]);
    vector<double> expectedvector2 = get_Linear_Transformation_expected_vector(dimension, plain_matrix, plain_matrix[0]);

    vector<double> errornorm1result(dimension);
    get_max_error_norm(plain_result, expectedvector2, dimension, errornorm1result);

}
