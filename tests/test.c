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
    setMaterialProps(&concrete, 100000, 38462);

    Section rectangle;
    setSectionProps(&rectangle, 0.01, 0.00001, 0.00001, 0.00001);

    Point p1;
    p1.x = 0;
    p1.y = 0;
    p1.z = 0;

    Point p2;
    p2.x = 0;
    p2.y = 0;
    p2.z = 2.44;

    Point v1;
    v1.x = 0;
    v1.y = 0.707;
    v1.z = 0.707;

    DegreesOfFreedom allFixed;
    setDegreesOfFreedomProps(&allFixed, false, false, false);

    DegreesOfFreedom allFree;
    setDegreesOfFreedomProps(&allFixed, true, true, true);

    NodalLoad load;

    Node n1;
    setNodeProps(&n1, p1, allFixed, allFixed, load);

    Node n2;
    setNodeProps(&n2, p2, allFree, allFree, load);

    FrameBar b1;
    setFrameBarProps(&b1, n1, n2, v1, &concrete, &rectangle);

    printMatrix(b1.stiffnessMatrix.global);

    return 0;
}
