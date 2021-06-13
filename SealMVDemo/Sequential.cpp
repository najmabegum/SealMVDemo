
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Common.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

void Sequential(int dimension, bool rescale)
{
    /*Set Seal context*/
    size_t poly_modulus_degree = 16384;
    /*size_t poly_modulus_degree = 8192;*/
    EncryptionParameters params(scheme_type::ckks);
    params.set_poly_modulus_degree(poly_modulus_degree);
    cout << "MAX BIT COUNT: " << CoeffModulus::MaxBitCount(poly_modulus_degree) << endl;    
    params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 40, 40, 40, 40, 60 }));
    /*params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 60 }));*/
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
    // --------------- MATRIX SET ------------------    
    cout << "Dimension Set 2: " << dimension << endl
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
            /*pod_matrix1_set2[i][j] = filler;
            pod_matrix2_set2[i][j] = static_cast<double>((j % 2) + 1);
            filler++;*/
        }
    }
    //double filler = 1.0;
    //int numberOfObs = dimension * dimension;
    //int maxMatrixEntry = numberOfObs / 2;
    //bool isMaxReached = false;
    //bool isReverseMaxEntryRepeat = false;
    //cout << "Matrix Maximum entry: " << maxMatrixEntry << endl
    //    << endl;

    //// Set 1
    //for (int i = 0; i < dimension; i++)
    //{
    //    for (int j = 0; j < dimension; j++)
    //    {
    //        isReverseMaxEntryRepeat = false;
    //        pod_matrix1_set2[i][j] = filler;
    //        pod_matrix2_set2[i][j] = static_cast<double>((j % 2) + 1);
    //        if (filler == maxMatrixEntry)
    //        {
    //            isMaxReached = true;
    //            isReverseMaxEntryRepeat = true;
    //        }
    //        if (isMaxReached)
    //        {
    //            filler--;
    //            if (filler < 0)
    //            {
    //                filler++;
    //            }
    //        }
    //        else
    //        {
    //            filler++;
    //        }
    //    }
    //}
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

    // Encrypt Diagonals    
    vector<Ciphertext> cipher_matrix1_set2(dimension), cipher_matrix2_set2(dimension);
    vector<Ciphertext> cipher_diagonal1_set2(dimension), cipher_diagonal2_set2(dimension);
    vector<Ciphertext> cipher_matrix1_lt_Seq(dimension);

    for (unsigned int i = 0; i < dimension; i++)
    {      
        encryptor.encrypt(plain_matrix1_set2[i], cipher_matrix1_set2[i]);
        encryptor.encrypt(plain_matrix2_set2[i], cipher_matrix2_set2[i]);
        encryptor.encrypt(plain_diagonal1_set2[i], cipher_diagonal1_set2[i]);
        encryptor.encrypt(plain_diagonal2_set2[i], cipher_diagonal2_set2[i]);
    }
    cout << "Encrypting Set 2 is Complete" << endl;   

    // ------------- SECOND COMPUTATION ----------------
    outf << "# index 1" << endl;
    outf << "# C_Vec . C_Mat" << endl;    

    /*Perform sequentially cipher matrix and cipher vector multiplication*/
    int lastIndex = dimension - 1;
   
    auto start_comp2_set2 = chrono::high_resolution_clock::now();
    //for (int x = 0; x < dimension; x++)
    //{
    //    if (x == 0)
    //    {
    //        encryptor.encrypt(plain_matrix1_set2[x], cipher_matrix1_set2[x]);
    //        cipher_matrix1_lt_Seq[x] = Linear_Transform_Cipher(cipher_matrix1_set2[x], cipher_diagonal1_set2, gal_keys, params, rescale);
    //    }
    //    else
    //    {
    //        /*Decrypt*/
    //        Plaintext pt_result_current;
    //        decryptor.decrypt(cipher_matrix1_lt_Seq[x-1], pt_result_current);

    //        // Decode
    //        vector<double> output_result_plain;
    //        ckks_encoder.decode(pt_result_current, output_result_plain);

    //        /*Comment two lines below for performance study*/
    //        cout << "Linear Transformation Set 2 Result:" << endl;
    //        print_partial_vector(output_result_plain, dimension);

    //        Ciphertext encodedCipherResult;
    //        Plaintext encodedPlainResult;
    //        ckks_encoder.encode(output_result_plain, scale, encodedPlainResult);
    //        encryptor.encrypt(encodedPlainResult, encodedCipherResult);
    //        if (cipher_diagonal1_set2[0].coeff_modulus_size() > encodedCipherResult.coeff_modulus_size())
    //        {
    //            for (unsigned int i = 0; i < dimension; i++)
    //            {
    //                evaluator.rescale_to_next_inplace(cipher_diagonal1_set2[i]);
    //            }
    //        }
    //        cipher_matrix1_lt_Seq[x] = Linear_Transform_Cipher(encodedCipherResult, cipher_diagonal1_set2, gal_keys, params, rescale);
    //    }
    //}

    encryptor.encrypt(plain_matrix1_set2[0], cipher_matrix1_set2[0]);
    Ciphertext finalresult = Linear_Transform_Cipher_Sequential_2(cipher_matrix1_set2[0], cipher_diagonal1_set2, gal_keys, params,
        rescale, sk, relin_keys,plain_diagonal1_set2, pk, all_diagonal1_set2, scale);

    //auto stop_comp2_set2 = chrono::high_resolution_clock::now();

    //auto duration_comp2_set2 = chrono::duration_cast<chrono::microseconds>(stop_comp2_set2 - start_comp2_set2);
    //cout << "\nTime to compute : " << duration_comp2_set2.count() << " microseconds" << endl;
    //outf << "100\t\t" << duration_comp2_set2.count() << endl;
   
    // /*Decrypt*/
    // int sizeArr = cipher_matrix1_lt_Seq.size();
    // Plaintext pt_result2_set2;
    // decryptor.decrypt(cipher_matrix1_lt_Seq[sizeArr-1], pt_result2_set2);

    // // Decode
    // vector<double> output_result2_set2;
    // ckks_encoder.decode(pt_result2_set2, output_result2_set2);

    // cout << "Linear Transformation Set 2 Result:" << endl;
    // print_partial_vector(output_result2_set2, dimension);

    // vector<double> expectedvector2 = get_Linear_Transformation_expected_sequentialVector(dimension, pod_matrix1_set2);
    // print_partial_vector(expectedvector2, dimension);
     /*get_max_error_norm_largeValues(output_result2_set2, expectedvector2, dimension);*/


    /*Decrypt*/
    
    Plaintext pt_result2_set2;
    decryptor.decrypt(finalresult, pt_result2_set2);

    // Decode
    vector<double> output_result2_set2;
    ckks_encoder.decode(pt_result2_set2, output_result2_set2);

    cout << "Linear Transformation Set 2 Result:" << endl;
    print_full_vector(output_result2_set2, dimension);

    /*vector<double> expectedvector2 = get_Linear_Transformation_expected_sequentialVector(dimension, pod_matrix1_set2);*/
    print_Linear_Transformation_expected_sequentialVector(dimension, pod_matrix1_set2);
    /*print_partial_vector(expectedvector2, dimension);*/
    /*get_max_error_norm(output_result2_set2, expectedvector2, dimension);*/

    outf.close();
}
