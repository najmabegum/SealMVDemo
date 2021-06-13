#include "Common.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

int main()
{
    cout << "Microsoft SEAL version: " << SEAL_VERSION << endl;
    while (true)
    {
        cout << "+---------------------------------------------------------+" << endl;
        cout << "| The following examples should be executed while reading |" << endl;
        cout << "| comments in associated files in native/examples/.       |" << endl;
        cout << "+---------------------------------------------------------+" << endl;
        cout << "| Examples                   |" << endl;
        cout << "+----------------------------+---------------------------------------------------------------------------------+" << endl;
        cout << "| 1.  Matriv Vector Multiplication - Dense (With rescale)                                                      |" << endl;
        cout << "| 2.  Matriv Vector Multiplication - Dense (Without rescale)                                                   |" << endl;
        cout << "| 3.  Matriv Vector Multiplication - Sparse (With rescale)                                                     |" << endl;
        cout << "| 4.  Matriv Vector Multiplication - Sparse (Without rescale)                                                  |" << endl;
        cout << "| 5.  Sequential Cipher Matriv Vector Multiplication - Dense (With rescale)                                    |" << endl;
        cout << "| 6.  Sequential Cipher Matriv Vector Multiplication - Dense (Without rescale)                                 |" << endl;
        cout << "| 7.  Matriv Vector Multiplication - Dense 1 matrix (Max Normalised and MSE) (With rescale)                    |" << endl;
        cout << "| 8.  Matriv Vector Multiplication - Dense 1 matrix (Max Normalised and MSE) (Without rescale)                 |" << endl;
        cout << "| 9.  Matriv Vector Multiplication - Dense (Max Normalised and MSE) (With rescale)                             |" << endl;
        cout << "| 10. Matriv Vector Multiplication - Dense (Max Normalised and MSE) (Without rescale)                          |" << endl;
        cout << "+----------------------------+---------------------------------------------------------------------------------+" << endl;
       
        size_t megabytes = MemoryManager::GetPool().alloc_byte_count() >> 20;
        cout << "[" << setw(10) << right << megabytes << " MB] "
            << "Total allocation from the memory pool" << endl;

        int selection = 0;
        bool valid = true;
        do
        {
            cout << endl << "> Run example (1 ~ 10) or exit (0): ";
            if (!(cin >> selection))
            {
                valid = false;
            }
            else if (selection < 0 || selection > 10)
            {
                valid = false;
            }
            else
            {
                valid = true;
            }
            if (!valid)
            {
                cout << "  [Beep~~] valid option: type 0 ~ 10" << endl;
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
            }
        } while (!valid);
        vector<int> dimensionValues = { 10,100,200,300,400,500,600,700,800,900,1000 };        
        switch (selection)
        {
        case 1:
            // input marix in increasing scale upto max value dimension*dimension
            MVMultiplicationDense(100,true);
            break;

        case 2:
            // input marix in increasing scale upto max value dimension*dimension
            MVMultiplicationDense(100, false);
            break;

        case 3:
            // input marix with main diagonal width 20
            MVMultiplicationSparse(100,0,19,true);
            break;

        case 4:
            // input marix with main diagonal width 20
            MVMultiplicationSparse(100, 0, 19, false);
            break;

        case 5:
            /*encrypted_result_1 = matrix * vector
            encrypted_result_2 = matrix * encrypted_result_1
            encrypted_result_3 = matrix * encrypted_result_2*/
            Sequential(7,true);
            break;        

        case 6:
            /*encrypted_result_1 = matrix * vector
            encrypted_result_2 = matrix * encrypted_result_1
            encrypted_result_3 = matrix * encrypted_result_2*/
            Sequential(7,false);
            break;

        case 7:
            // input marix is a ones matrix            
            for (int i = 0;i < dimensionValues.size();i++)
            {
                cout << endl << "> Current Dimension: " << dimensionValues[i] << endl;
                MVDenseMatrix(dimensionValues[i], true);
            }
            break;

        case 8:
            // input marix is a ones matrix
            for (int i = 0;i < dimensionValues.size();i++)
            {
                cout << endl << "> Current Dimension: " << dimensionValues[i] << endl;
                MVDenseMatrix(dimensionValues[i], false);
            }                      
            break;

        case 9:            
            // input marix increases by 1 until dimension * dimension /2, then decreses by 1
            for (int i = 0;i < dimensionValues.size();i++)
            {
                cout << endl << "> Current Dimension: " << dimensionValues[i] << endl;
                Mult_MVDense(dimensionValues[i], true);
            }            
            break;

        case 10:
            // input marix increases by 1 until dimension * dimension /2, then decreses by 1
            for (int i = 0;i < dimensionValues.size();i++)
            {
                cout << endl << "> Current Dimension: " << dimensionValues[i] << endl;
                Mult_MVDense(dimensionValues[i], false);
            }            
            break;

        case 0:
            return 0;
        }
    }

    return 0;
}
