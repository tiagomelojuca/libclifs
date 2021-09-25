#include <stdio.h>
#include <stdlib.h>

#include "../src/libclifs.h"

#define A  0.01
#define I  0.00001
#define E  100000000.0
#define G  38462000.0

#define VX 0.000
#define VY 0.707
#define VZ 0.707

int main()
{
    Point    avec = createPoint(VX, VY, VZ);
    Section  rect = createSection(A, I, I, I);
    Material conc = createMaterial(E, G);

    DegreesOfFreedom allFree = createDegreesOfFreedom(true, true, true);
    DegreesOfFreedom allFix = createDegreesOfFreedom(false, false, false);

    NodalLoad l1 = createNodalLoad(0.000,  0.000,   0.000,  0.000,  0.000,  0.000);
    NodalLoad l2 = createNodalLoad(8.900,  0.000,   0.000,  0.000,  0.000,  0.000);
    NodalLoad l3 = createNodalLoad(0.000, 17.800,   0.000,  0.000,  0.000,  0.000);
    NodalLoad l4 = createNodalLoad(0.000,  0.000,  -4.450,  0.000,-10.858,  0.000);
    NodalLoad l5 = createNodalLoad(0.000,  0.000,   0.000,  0.000,  0.000,  0.000);

    Node n1 = createNode(createPoint(0.000,  0.000,  0.000), allFix,  allFix,  l1);
    Node n2 = createNode(createPoint(0.000,  0.000,  2.440), allFree, allFree, l2);
    Node n3 = createNode(createPoint(2.440,  0.000,  2.440), allFree, allFree, l3);
    Node n4 = createNode(createPoint(4.880,  0.000,  2.440), allFree, allFree, l4);
    Node n5 = createNode(createPoint(7.320,  2.440,  0.000), allFix,  allFix,  l5);

    GlobalSystem* g = createGlobalSystem();
    
    insertNodeGlobalSystem(g, n1);
    insertNodeGlobalSystem(g, n2);
    insertNodeGlobalSystem(g, n3);
    insertNodeGlobalSystem(g, n4);
    insertNodeGlobalSystem(g, n5);

    FrameBar b1 = createFrameBar(getNode(g, 1), getNode(g, 2), avec, &conc, &rect);
    FrameBar b2 = createFrameBar(getNode(g, 2), getNode(g, 3), avec, &conc, &rect);
    FrameBar b3 = createFrameBar(getNode(g, 3), getNode(g, 4), avec, &conc, &rect);
    FrameBar b4 = createFrameBar(getNode(g, 4), getNode(g, 5), avec, &conc, &rect);

    insertFrameBarGlobalSystem(g, b1);
    insertFrameBarGlobalSystem(g, b2);
    insertFrameBarGlobalSystem(g, b3);
    insertFrameBarGlobalSystem(g, b4);

    mountGlobalSystem(g);

    printf("SUPPORT REACTIONS:\n");
    for(int i = 0; i < g->numEqConstraints; i++) {
        printf("%.2f\n", g->vecSupportReactions[i][0]);
    }

    freeGlobalSystem(g);
    return 0;
}
