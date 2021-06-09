
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

	// --------------- MATRIX SET ------------------

	cout << "Dimension of Matrix : " << dimension << endl
		<< endl;

	vector<vector<double>> pod_matrix1_set2(dimension, vector<double>(dimension));
	vector<vector<double>> pod_matrix2_set2(dimension, vector<double>(dimension));

	// Fill input matrices
	double filler = 1.0;

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
	// test print output
	/*cout << "Diagonal Set 2" << endl;
	cout << "\nDiagonal Set 2 Result:" << endl;
	for (unsigned int i = 0; i < dimension; i++)
	{
		cout << "\t[";
		for (unsigned int j = 0; j < dimension-1; j++)
		{
			cout << all_diagonal1_set2[i][j] << ", ";
		}
		cout << all_diagonal1_set2[i][dimension-1];
		cout << "]" << endl;
	}
	cout << "\n"
		<< endl;*/


		// Test LinearTransform here    

			///*UNCOMMENT TO DEBUG LINEAR TRANSFORM FUNCTION*/
			//// Fill ct
			//Ciphertext ct_rotated;
			//evaluator.rotate_vector(cipher_matrix1_set2[0], dimension, gal_keys, ct_rotated);
			//Ciphertext ct;
			//evaluator.add(cipher_matrix1_set2[0], ct_rotated, ct);
			//// Add epsilon to avoid negative numbers
			//vector<double> epsilon_vec(poly_modulus_degree / 2);
			//for (int i = 0; i < epsilon_vec.size(); i++)
			//{
			//    epsilon_vec[i] = 0.0000;
			//}
			//Plaintext epsilon_plain;
			//ckks_encoder.encode(epsilon_vec, scale, epsilon_plain);
			//evaluator.add_plain_inplace(ct, epsilon_plain);
			//// test fill ct
			//Plaintext test_fill;
			//decryptor.decrypt(ct, test_fill);
			//vector<double> out_fill;
			//ckks_encoder.decode(test_fill, out_fill);
			//cout << "Filled CT:\n"
			//    << endl;
			//for (int i = 0; i < dimension; i++)
			//{
			//    cout << out_fill[i] << ", ";
			//}
			//cout << "\n"
			//    << endl;
			//Ciphertext ct_prime;
			//// ct` = CMult(ct, u0)
			//evaluator.multiply_plain(ct, plain_diagonal1_set2[0], ct_prime);
			//// test mult plain 0
			//Plaintext test_0;
			//decryptor.decrypt(ct_prime, test_0);
			//vector<double> out_test_0;
			//ckks_encoder.decode(test_0, out_test_0);
			//cout << "CT_Prime 0 :\n"
			//    << endl;
			//for (int i = 0; i < dimension; i++)
			//{
			//    cout << out_test_0[i] << ", ";
			//}
			//cout << "\n"
			//    << endl;
			//for (int l = 1; l < dimension; l++)
			//{
			//    // ct` = Add(ct`, CMult(Rot(ct, l), ul))
			//    Ciphertext temp_rot;
			//    Ciphertext temp_mul;
			//    evaluator.rotate_vector(ct, l, gal_keys, temp_rot);
			//    evaluator.multiply_plain(temp_rot, plain_diagonal1_set2[l], temp_mul);
			//    evaluator.add_inplace(ct_prime, temp_mul);
			//    // test decrypt
			//    Plaintext temp_rot_plain;
			//    Plaintext temp_mul_plain;
			//    Plaintext temp_ct_prime;
			//    decryptor.decrypt(temp_rot, temp_rot_plain);
			//    decryptor.decrypt(temp_mul, temp_mul_plain);
			//    decryptor.decrypt(ct_prime, temp_ct_prime);
			//    // test decode
			//    vector<double> test_out_rot, test_out_mul, test_ct_prime;
			//    vector<double> test_diag;
			//    ckks_encoder.decode(temp_ct_prime, test_ct_prime);
			//    ckks_encoder.decode(temp_mul_plain, test_out_mul);
			//    ckks_encoder.decode(temp_rot_plain, test_out_rot);
			//    ckks_encoder.decode(plain_diagonal1_set2[l], test_diag);
			//    cout << "Rotation " << l << "\n"
			//        << endl;
			//    cout << "\nrotated vec:\n\t[";
			//    for (int j = 0; j < dimension; j++)
			//    {
			//        cout << test_out_rot[j] << ", ";
			//    }
			//    cout << "\nDiagonal vec:\n\t[";
			//    for (int j = 0; j < dimension; j++)
			//    {
			//        cout << test_diag[j] << ", ";
			//    }
			//    cout << "\nMult vec vec:\n\t[";
			//    for (int j = 0; j < dimension; j++)
			//    {
			//        cout << test_out_mul[j] << ", ";
			//    }
			//    cout << "\nCt_prime vec:\n\t[";
			//    for (int j = 0; j < dimension; j++)
			//    {
			//        cout << test_ct_prime[j] << ", ";
			//    }
			//    cout << "\n"
			//        << endl;
			//}


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
	vector<double> errornorm1result(dimension);
	get_max_error_norm(output_result_set2, expectedvector1, dimension, errornorm1result);

	// ------------- SECOND COMPUTATION ----------------
	outf << "# index 1" << endl;
	outf << "# C_Vec . C_Mat" << endl;


	/*Perform Linear Transform on Cipher*/

	auto start_comp2_set2 = chrono::high_resolution_clock::now();
	Ciphertext ct_prime2_set2 = Linear_Transform_Cipher(cipher_matrix1_set2[0], cipher_diagonal1_set2, gal_keys, params, rescale);
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
	vector<double> errornorm2result(dimension);
	get_max_error_norm(output_result2_set2, expectedvector2, dimension, errornorm2result);

	outf.close();
}