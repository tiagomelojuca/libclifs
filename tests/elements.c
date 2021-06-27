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

    FrameBar b1 = createFrameBar(
        createNode(createPoint(0.0, 0.0, 0.0),
                   createDegreesOfFreedom(false, false, false),
                   createDegreesOfFreedom(false, false, false),
                   createNodalLoad(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        createNode(createPoint(0.0, 0.0, 2.44),
                   createDegreesOfFreedom(true, true, true),
                   createDegreesOfFreedom(true, true, true),
                   createNodalLoad(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        createPoint(0.0, 0.707, 0.707),
        &concrete,
        &rectangle
    );

    printMatrix(b1.stiffnessMatrix.global);

    return 0;
}
