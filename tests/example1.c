#include <stdio.h>
#include <stdlib.h>
#include "../src/libclifs.h"

void printMatrix(double** matrix, int nRows, int nColumns)
{
    for(int i = 0; i < nRows; i++) {
        for(int j = 0; j < nColumns; j++) {
            printf("%.1f  ", matrix[i][j]);
        } printf("\n");
    }
}

int main()
{
    GlobalSystem g;
    initGlobalSystem(&g);

    Material concrete = createMaterial(100000.0, 38462.0);
    Section rectangle = createSection(0.01, 0.00001, 0.00001, 0.00001);

    DegreesOfFreedom allFree = createDegreesOfFreedom(true, true, true);
    DegreesOfFreedom allFix = createDegreesOfFreedom(false, false, false);

    NodalLoad l1 = createNodalLoad(0.0,  0.0,   0.0, 0.0,     0.0, 0.0);
    NodalLoad l2 = createNodalLoad(8.9,  0.0,   0.0, 0.0,     0.0, 0.0);
    NodalLoad l3 = createNodalLoad(0.0, 17.8,   0.0, 0.0,     0.0, 0.0);
    NodalLoad l4 = createNodalLoad(0.0,  0.0, -4.45, 0.0, -10.858, 0.0);
    NodalLoad l5 = createNodalLoad(0.0,  0.0,   0.0, 0.0,     0.0, 0.0);

    Node n1 = createNode(createPoint( 0.0,  0.0,  0.0), allFix,  allFix,  l1);
    Node n2 = createNode(createPoint( 0.0,  0.0, 2.44), allFree, allFree, l2);
    Node n3 = createNode(createPoint(2.44,  0.0, 2.44), allFree, allFree, l3);
    Node n4 = createNode(createPoint(4.88,  0.0, 2.44), allFree, allFree, l4);
    Node n5 = createNode(createPoint(7.32, 2.44,  9.0), allFix,  allFix,  l5);
    
    insertNodeGlobalSystem(&g, n1);
    insertNodeGlobalSystem(&g, n2);
    insertNodeGlobalSystem(&g, n3);
    insertNodeGlobalSystem(&g, n4);
    insertNodeGlobalSystem(&g, n5);

    Point auxvec = createPoint(0.0, 0.707, 0.707);
    FrameBar b1 = createFrameBar(n1, n2, auxvec, &concrete, &rectangle);
    FrameBar b2 = createFrameBar(n2, n3, auxvec, &concrete, &rectangle);
    FrameBar b3 = createFrameBar(n3, n4, auxvec, &concrete, &rectangle);
    FrameBar b4 = createFrameBar(n4, n5, auxvec, &concrete, &rectangle);

    insertFrameBarGlobalSystem(&g, b1);
    insertFrameBarGlobalSystem(&g, b2);
    insertFrameBarGlobalSystem(&g, b3);
    insertFrameBarGlobalSystem(&g, b4);

    printf("Number of Nodes: %i\n", g.nodeArray.used);
    printf("Number of FrameBars: %i\n", g.framebarsArray.used);
    printf("Number of Equations: %i\n", g.numEquations);
    printf("Number of Freedoms: %i\n", g.numEqFreedoms);
    printf("Number of Constraints: %i\n", g.numEqConstraint);

    mountGlobalSystem(&g);
    printMatrix(g.mtxStiffness, g.numEquations, g.numEquations); printf("\n");
    printMatrix(g.mtxNodalLoads, g.numEquations, 1); printf("\n");
    printMatrix(g.mtxDisplacements, g.numEquations, 1); printf("\n");
    printMatrix(g.mtxFreedoms, DOG, g.nodeArray.used); printf("\n");
    printMatrix(g.mtxSpreading, g.framebarsArray.used, SM); printf("\n");
    
    freeGlobalSystem(&g);
    return 0;
}
