#include "Common.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

int main()
{
    cout << "Microsoft SEAL version: " << SEAL_VERSION << endl;


    cout << "Initiating the Sequential Procedure " << endl;

    SequentialMV(6, 100, true);

    return 0;
}
