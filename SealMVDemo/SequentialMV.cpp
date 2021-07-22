
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Common.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

void SequentialMV(int iterations, int dimension,bool rescale)
{
    /*Set Seal context*/
    EncryptionParameters params(scheme_type::ckks);

    /*size_t poly_modulus_degree = 8192;*/
    size_t poly_modulus_degree = 16384;
    params.set_poly_modulus_degree(poly_modulus_degree);
    cout << "MAX BIT COUNT: " << CoeffModulus::MaxBitCount(poly_modulus_degree) << endl;
    
    /*params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 60 }));*/   
    params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 60 }));
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
    double scaleBit = 30;
    double scale = pow(2.0, 30);

    
  
    cout << "Dimension of the matrix: " << dimension << endl
        << endl;

    vector<vector<double>> plain_matrix(dimension, vector<double>(dimension));    

    // Filling the matrix with 1.0, 2.0 and so on... And in the middle we start decrementing the filler again...
    double filler = 1.0;

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            plain_matrix[i][j] = filler;  

            if (i >= (int) dimension / 2)
            {
                filler--;
            }
            else
            {
                filler++;
            }
        }
    }
    
    cout << "This is the matrix that will be used in the multiplication: " << endl;

    print_partial_matrix(plain_matrix);    

    // Get all diagonals    
    vector<vector<double>> matrix_diagonals(dimension, vector<double>(dimension));    

    for (int i = 0; i < dimension; i++)
    {
        matrix_diagonals[i] = get_diagonal(i, plain_matrix);        
    }

    cout << "These are the diagonals of the matrix:" << endl;
    print_partial_matrix(matrix_diagonals);   

    // Encode Matrices into vectors with Diagonals   

    // Encoding the matrix and the diagonals, should we actually encode the matrix? Is it needed? Don't we lose time with it??
    // Yes we do lose time with it, since at the end we only get one row of this encoded/encrypted matrix.
    // So what do we do?
    // We encode and encrypt ONLY the first row that we use!

    Plaintext encoded_vector;
    Ciphertext encrypted_vector;
    // We get the first row of the matrix and use it as a vector for the multiplication 

    cout << "Encoding and encrypting the vector..." << endl;

    ckks_encoder.encode(plain_matrix[0], scale, encoded_vector);
    encryptor.encrypt(encoded_vector, encrypted_vector);

    cout << "Encryption of the vector complete." << endl;


    vector<Plaintext> encoded_matrix_diagonals(dimension);

    cout << "Encoding Matrix (diagonals)...." << endl;

    for (int i = 0; i < dimension; i++)
    {
        ckks_encoder.encode(matrix_diagonals[i], scale, encoded_matrix_diagonals[i]);        
    }

    cout << "Encoding of the matrix (diagonals) Complete." << endl;

    // Encrypt Diagonals    
    vector<Ciphertext> encrypted_matrix_diagonals(dimension);
    vector<Ciphertext> cipher_matrix1_lt_Seq(dimension);

    cout << "Encrypting Matrix (diagonals)...." << endl;

    for (unsigned int i = 0; i < dimension; i++)
    {
        encryptor.encrypt(encoded_matrix_diagonals[i], encrypted_matrix_diagonals[i]);        
    }
    cout << "Encrypting Matrix (diagonals) complete." << endl;


    // ---------------- THE PART THAT IS DIFFERENT FROM THE OTHER IMPLEMENTATIONS ----------------

    cout << "Starting to perform sequential operations." << endl;


    int lastIndex = dimension - 1;

    auto start_of_sequential_ops = chrono::high_resolution_clock::now();    

    Ciphertext finalresult = Linear_Transform_Cipher_Sequential_Dense(encrypted_vector, encrypted_matrix_diagonals, gal_keys, params,
        rescale, sk, relin_keys, iterations, context, scaleBit);

    auto end_of_sequential_ops = chrono::high_resolution_clock::now();

    auto duration_of_sequential_ops = chrono::duration_cast<chrono::microseconds>(end_of_sequential_ops - start_of_sequential_ops);

    cout << "\nTotal Time to compute all the sequential operations : " << duration_of_sequential_ops.count()/ 1000000 << " Seconds" << endl;

    print_Linear_Transformation_expected_sequential_dense(iterations,dimension, plain_matrix);
}
