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
    NodeArray a;
    initNodeArray(&a, 1);

    insertNodeArray(&a, createNode(createPoint(0.0, 0.0, 0.0),
                                   createDegreesOfFreedom(false, false, false),
                                   createDegreesOfFreedom(false, false, false),
                                   createNodalLoad(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));

    insertNodeArray(&a, createNode(createPoint(0.0, 0.0, 2.44),
                                   createDegreesOfFreedom(true, true, true),
                                   createDegreesOfFreedom(true, true, true),
                                   createNodalLoad(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));

    printf("%f", a.nodes[1].position.z);

    return 0;
}
