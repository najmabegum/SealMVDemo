
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Common.h"
#include "seal/seal.h"


using namespace std;
using namespace seal;

void MVDenseMatrix(int dimension, bool rescale)
{
	/*Setting up Seal Context*/
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
	// cout << "Coeff Modulus Back Value: " << params.coeff_modulus().back().value() << endl;
	double scale = pow(2.0, 40);


	// --------------- Create vectors to hold normalised error values------------------

	int numberOfObs = dimension * dimension;
	vector<double> error_result(dimension);
	vector<double> error_maxNormalised_result(dimension);
	vector<double> error_mse_result(dimension);


	vector<vector<double>> plain_matrix(dimension, vector<double>(dimension));

	// Fill input matrices
	double filler = 1.0;
	int maxMatrixEntry = numberOfObs / 2;
	bool isMaxReached = false;
	bool isReverseMaxEntryRepeat = false;

	// cout << "Matrix Maximum entry: " << maxMatrixEntry << endl << endl;

	// Set 1
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			plain_matrix[i][j] = filler;
		}
	}

	// print_partial_matrix(plain_matrix);


	// Get all diagonals        
	vector<vector<double>> plain_diagonals(dimension, vector<double>(dimension));

	for (int i = 0; i < dimension; i++)
	{
		plain_diagonals[i] = get_diagonal(i, plain_matrix);
	}

	/*cout << "Diagonals of the Matrix: " << endl;
	print_partial_matrix(plain_diagonals);*/


	// Encode Matrices into vectors with Diagonals        
	vector<Plaintext> encoded_diagonals(dimension);

	for (int i = 0; i < dimension; i++)
	{
		ckks_encoder.encode(plain_diagonals[i], scale, encoded_diagonals[i]);
	}

	//cout << "Encoding Complete" << endl;

	// Encrypt the matrices with Diagonals    
	vector<Ciphertext> encrypted_diagonals(dimension);

	for (unsigned int i = 0; i < dimension; i++)
	{
		encryptor.encrypt(encoded_diagonals[i], encrypted_diagonals[i]);
	}
	//cout << "Encrypting Complete" << endl;


	//cout << "Encoding and encrypting the Vector to be used for the multiplication" << endl;

	Plaintext encoded_vector;
	Ciphertext encrypted_vector;
	ckks_encoder.encode(plain_matrix[0], scale, encoded_vector);
	encryptor.encrypt(encoded_vector, encrypted_vector);

	//cout << "Encoding and Encryption of the Vector complete.";

	/*Perform Linear Transform on Cipher*/

	auto start_comp2_set2 = chrono::high_resolution_clock::now();
	Ciphertext ct_prime2_set2 = Linear_Transform_Cipher(encrypted_vector, encrypted_diagonals, gal_keys, params, rescale, relin_keys);
	auto stop_comp2_set2 = chrono::high_resolution_clock::now();

	// Decrypt
	Plaintext encoded_results;
	decryptor.decrypt(ct_prime2_set2, encoded_results);

	// Decode
	vector<double> plain_results;
	ckks_encoder.decode(encoded_results, plain_results);

	auto duration_comp2_set2 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set2 - start_comp2_set2);

	cout << "HE Result:" << endl;
	print_partial_vector(plain_results, dimension);

	// Check result
	cout << "Expected Result: " << endl;
	test_Linear_Transformation(dimension, plain_matrix, plain_matrix[0]);

	vector<double> expectedvector2 = get_Linear_Transformation_expected_vector(dimension, plain_matrix, plain_matrix[0]);

	cout << "\nTime to compute C_vec . C_mat: " << endl << 
		duration_comp2_set2.count() << " microseconds" << endl;

	// Calculate Error Norms
	error_maxNormalised_result[1] = get_max_error_norm_normalised(maxMatrixEntry, plain_results, expectedvector2, dimension, error_result);
	double mseSummationValue2 = get_mse(plain_results, expectedvector2, dimension, error_result);
	error_mse_result[1] = mseSummationValue2 / (dimension * dimension);

	/*Print Error Norms*/

	cout << "Max Normalised Error: " << endl <<
		error_maxNormalised_result[1] << endl;
	
	cout << "MSE Norm: " << endl <<
		error_mse_result[1] << endl;

}