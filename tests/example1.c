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

    Frame frame = createFrame();

    pushFrameNode(frame, createPoint(0.000,  0.000,  0.000), allFix,  allFix,  l1);
    pushFrameNode(frame, createPoint(0.000,  0.000,  2.440), allFree, allFree, l2);
    pushFrameNode(frame, createPoint(2.440,  0.000,  2.440), allFree, allFree, l3);
    pushFrameNode(frame, createPoint(4.880,  0.000,  2.440), allFree, allFree, l4);
    pushFrameNode(frame, createPoint(7.320,  2.440,  0.000), allFix,  allFix,  l5);

    pushFrameBar(frame, 1, 2, avec, &conc, &rect);
    pushFrameBar(frame, 2, 3, avec, &conc, &rect);
    pushFrameBar(frame, 3, 4, avec, &conc, &rect);
    pushFrameBar(frame, 4, 5, avec, &conc, &rect);

    solveFrame(frame);

    printf("SUPPORT REACTIONS:\n");
    for(int i = 0; i < frame->numEqConstraints; i++) {
        printf("%.2f\n", frame->vecSupportReactions[i][0]);
    }

    freeFrame(frame);
    return 0;
}
