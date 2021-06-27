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
    Material concrete = createMaterial(100000.0, 38462.0);
    Section rectangle = createSection(0.01, 0.00001, 0.00001, 0.00001);

    Point p1 = createPoint(0.0,   0.0,   0.0);
    Point p2 = createPoint(0.0,   0.0,  2.44);
    Point v1 = createPoint(0.0, 0.707, 0.707);

    DegreesOfFreedom allFixed = createDegreesOfFreedom(false, false, false);
    DegreesOfFreedom allFree = createDegreesOfFreedom(true, true, true);;

    NodalLoad load = createNodalLoad(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    Node n1 = createNode(p1, allFixed, allFixed, load);
    Node n2 = createNode(p2, allFree, allFree, load);

    FrameBar b1 = createFrameBar(n1, n2, v1, &concrete, &rectangle);

    printMatrix(b1.stiffnessMatrix.global);

    return 0;
}
