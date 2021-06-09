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
        cout << "+----------------------------+-----------------------------------------------------+" << endl;
        cout << "| 1. Matriv Vector Multiplication - Dense (With rescale)                           |" << endl;
        cout << "| 2. Matriv Vector Multiplication - Dense (Without rescale)                        |" << endl;
        cout << "| 3. Matriv Vector Multiplication - Sparse (With rescale)                          |" << endl;
        cout << "| 4. Matriv Vector Multiplication - Sparse (Without rescale)                       |" << endl;
        cout << "| 5. Sequential Cipher Matriv Vector Multiplication - Dense (With rescale)         |" << endl;
        cout << "| 6. Sequential Cipher Matriv Vector Multiplication - Dense (Without rescale)      |" << endl;
        cout << "| 7. Matriv Vector Multiplication - Dense 1 matrix (With rescale)                           |" << endl;
        cout << "| 8. Matriv Vector Multiplication - Dense 1 matrix (Without rescale)                           |" << endl;
        cout << "+----------------------------+----------------------------+" << endl;
       
        size_t megabytes = MemoryManager::GetPool().alloc_byte_count() >> 20;
        cout << "[" << setw(6) << right << megabytes << " MB] "
            << "Total allocation from the memory pool" << endl;

        int selection = 0;
        bool valid = true;
        do
        {
            cout << endl << "> Run example (1 ~ 6) or exit (0): ";
            if (!(cin >> selection))
            {
                valid = false;
            }
            else if (selection < 0 || selection > 6)
            {
                valid = false;
            }
            else
            {
                valid = true;
            }
            if (!valid)
            {
                cout << "  [Beep~~] valid option: type 0 ~ 6" << endl;
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
            }
        } while (!valid);

        switch (selection)
        {
        case 1:
            MVMultiplicationDense(100,true);
            break;

        case 2:
            MVMultiplicationDense(100, false);
            break;

        case 3:
            MVMultiplicationSparse(100,0,19,true);
            break;

        case 4:
            MVMultiplicationSparse(100, 0, 19, false);
            break;

        case 5:
            Sequential(10,true);
            break;        

        case 6:
            Sequential(10,false);
            break;

        case 7:
            MVDenseMatrix(10, true);
            break;

        case 8:
            MVDenseMatrix(10, false);
            break;

        case 0:
            return 0;
        }
    }

    return 0;
}
