#include <iostream>
#include "matrix.h"


int main()
{
    float data1[]{1, 2, 3, 1.03, 9, 2, 5, 3, 1, 0.2, 1, 9.1, 4.89, 1.23, 5, 1};
    matrix A(4, 4, data1);
    A.print();

    matrix B = !A;
    B.print();

    A.eye();
    A.print();

    return 0;
}
