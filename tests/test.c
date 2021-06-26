#include <stdio.h>
#include "../src/libclifs.h"

void printMatrix(double _matrix[12][12]) {
    for(int i = 0; i < 12; i++) {
        for(int j = 0; j < 12; j++) {
            printf("%f    ", _matrix[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    Material concrete;
    concrete.e = 100000;
    concrete.g = 38462;

    Section rectangle;
    rectangle.a = 0.01;
    rectangle.iy = 0.00001;
    rectangle.iz = 0.00001;
    rectangle.j = 0.00001;

    Node n1;
    n1.position.x = 0;
    n1.position.y = 0;
    n1.position.z = 0;

    Node n2;
    n2.position.x = 0;
    n2.position.y = 0;
    n2.position.z = 2.44;

    Point v1;
    v1.x = 0;
    v1.y = 0.707;
    v1.z = 0.707;

    Bar b1;
    setBarProps(&b1, n1, n2, v1, &concrete, &rectangle);

    StiffnessMatrix m1;
    setStiffnessMatrix(&m1, &b1);

    printMatrix(m1.global);

    return 0;
}
