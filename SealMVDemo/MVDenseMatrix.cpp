
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

	// --------------- Create vectors to hold normalised error values------------------

	int numberOfObs = dimension * dimension;
	vector<double> error_result(dimension);
	vector<double> error_maxNormalised_result(dimension);
	vector<double> error_mse_result(dimension);

	// --------------- MATRIX SET ------------------

	cout << "Dimension of Matrix : " << dimension << endl
		<< endl;

	vector<vector<double>> pod_matrix1_set2(dimension, vector<double>(dimension));
	vector<vector<double>> pod_matrix2_set2(dimension, vector<double>(dimension));

	// Fill input matrices
	double filler = 1.0;
	int maxMatrixEntry = numberOfObs / 2;
	bool isMaxReached = false;
	bool isReverseMaxEntryRepeat = false;
	cout << "Matrix Maximum entry: " << maxMatrixEntry << endl
		<< endl;

	// Set 1
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			pod_matrix1_set2[i][j] = filler;
			pod_matrix2_set2[i][j] = filler;
			/*pod_matrix2_set2[i][j] = static_cast<double>((j % 2) + 1);*/
			//filler++;
		}
	}
	print_partial_matrix(pod_matrix1_set2);
	print_partial_matrix(pod_matrix2_set2);


	// Get all diagonals        
	vector<vector<double>> all_diagonal1_set2(dimension, vector<double>(dimension));
	vector<vector<double>> all_diagonal2_set2(dimension, vector<double>(dimension));

	for (int i = 0; i < dimension; i++)
	{
		all_diagonal1_set2[i] = get_diagonal(i, pod_matrix1_set2);
		all_diagonal2_set2[i] = get_diagonal(i, pod_matrix2_set2);
	}

	cout << "Diagonal 1 Set 2 Expected:" << endl;
	print_partial_matrix(all_diagonal1_set2);

	cout << "Diagonal 2 Set 2 Expected:" << endl;
	print_partial_matrix(all_diagonal2_set2);

	// Encode Matrices into vectors with Diagonals        
	vector<Plaintext> plain_matrix1_set2(dimension), plain_matrix2_set2(dimension);
	vector<Plaintext> plain_diagonal1_set2(dimension), plain_diagonal2_set2(dimension);

	for (int i = 0; i < dimension; i++)
	{
		ckks_encoder.encode(pod_matrix1_set2[i], scale, plain_matrix1_set2[i]);
		ckks_encoder.encode(pod_matrix2_set2[i], scale, plain_matrix2_set2[i]);
		ckks_encoder.encode(all_diagonal1_set2[i], scale, plain_diagonal1_set2[i]);
		ckks_encoder.encode(all_diagonal2_set2[i], scale, plain_diagonal2_set2[i]);
	}

	cout << "Encoding Set 2 is Complete" << endl;

	// Encrypt the matrices with Diagonals    
	vector<Ciphertext> cipher_matrix1_set2(dimension), cipher_matrix2_set2(dimension);
	vector<Ciphertext> cipher_diagonal1_set2(dimension), cipher_diagonal2_set2(dimension);

	for (unsigned int i = 0; i < dimension; i++)
	{
		encryptor.encrypt(plain_matrix1_set2[i], cipher_matrix1_set2[i]);

		encryptor.encrypt(plain_matrix2_set2[i], cipher_matrix2_set2[i]);

		encryptor.encrypt(plain_diagonal1_set2[i], cipher_diagonal1_set2[i]);

		encryptor.encrypt(plain_diagonal2_set2[i], cipher_diagonal2_set2[i]);

	}
	cout << "Encrypting Set 2 is Complete" << endl;

	// ------------- FIRST COMPUTATION ----------------
	outf << "# index 0" << endl;
	outf << "# C_Vec . P_Mat" << endl;

	// Testing Decryption
	cout << "Decrypting Set 2" << endl;
	for (unsigned int i = 0; i < dimension; i++)
	{
		decryptor.decrypt(cipher_matrix1_set2[i], plain_matrix1_set2[i]);
		decryptor.decrypt(cipher_matrix2_set2[i], plain_matrix2_set2[i]);
		decryptor.decrypt(cipher_diagonal1_set2[i], plain_diagonal1_set2[i]);
		decryptor.decrypt(cipher_diagonal2_set2[i], plain_diagonal2_set2[i]);
	}
	// test decode here
	// test decrypt here

	// Testing Decoding
	cout << "Decoding Set 2" << endl;
	for (unsigned int i = 0; i < dimension; i++)
	{
		ckks_encoder.decode(plain_diagonal1_set2[i], all_diagonal1_set2[i]);
		ckks_encoder.decode(plain_diagonal2_set2[i], all_diagonal2_set2[i]);
	}

	/*Perform Linear Transform on Plain*/

	auto start_comp1_set2 = chrono::high_resolution_clock::now();
	Ciphertext ct_prime_set2 = Linear_Transform_Plain(cipher_matrix1_set2[0], plain_diagonal1_set2, gal_keys, params, rescale);
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
	print_partial_vector(output_result_set2, dimension);

	// Check result
	cout << "Expected output Set 2: " << endl;
	test_Linear_Transformation(dimension, pod_matrix1_set2, pod_matrix1_set2[0]);
	vector<double> expectedvector1 = get_Linear_Transformation_expected_vector(dimension, pod_matrix1_set2, pod_matrix1_set2[0]);
	
	// Calculate Error Norms
	error_maxNormalised_result[0] = get_max_error_norm_normalised(maxMatrixEntry, output_result_set2, expectedvector1, dimension, error_result);
	double mseSummationValue1 = get_mse(output_result_set2, expectedvector1, dimension, error_result);
	error_mse_result[0] = mseSummationValue1 / (dimension * dimension);

	// ------------- SECOND COMPUTATION ----------------
	outf << "# index 1" << endl;
	outf << "# C_Vec . C_Mat" << endl;


	/*Perform Linear Transform on Cipher*/

	auto start_comp2_set2 = chrono::high_resolution_clock::now();
	Ciphertext ct_prime2_set2 = Linear_Transform_Cipher(cipher_matrix1_set2[0], cipher_diagonal1_set2, gal_keys, params, rescale,relin_keys);
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
	print_partial_vector(output_result2_set2, dimension);

	// Check result
	cout << "Expected output Set 2: " << endl;
	test_Linear_Transformation(dimension, pod_matrix1_set2, pod_matrix1_set2[0]);
	vector<double> expectedvector2 = get_Linear_Transformation_expected_vector(dimension, pod_matrix1_set2, pod_matrix1_set2[0]);
	// Calculate Error Norms
	error_maxNormalised_result[1] = get_max_error_norm_normalised(maxMatrixEntry, output_result2_set2, expectedvector2, dimension, error_result);
	double mseSummationValue2 = get_mse(output_result2_set2, expectedvector2, dimension, error_result);
	error_mse_result[1] = mseSummationValue1 / (dimension * dimension);

	/*Print Error Norms*/
	cout << "Max Normalised Error:" << endl;
	cout << "C_Vec . P_Mat:" << error_maxNormalised_result[0] << endl;
	cout << "C_Vec . C_Mat:" << error_maxNormalised_result[1] << endl;

	cout << "MSE Norm:" << endl;
	cout << "C_Vec . P_Mat:" << error_mse_result[0] << endl;
	cout << "C_Vec . C_Mat:" << error_mse_result[1] << endl;

	outf.close();
}