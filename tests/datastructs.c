#include <stdio.h>
#include "../src/libclifs.h"

int main()
{
    Material concrete = createMaterial(100000.0, 38462.0);
    Section rectangle = createSection(0.01, 0.00001, 0.00001, 0.00001);

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

    FrameBarArray fb;
    initFrameBarArray(&fb, 1);

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

    FrameBar b2 = createFrameBar(
        createNode(createPoint(1.0, 2.0, 0.0),
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

    insertFrameBarArray(&fb, b1);
    insertFrameBarArray(&fb, b2);

    printf("%f", fb.framebars[1].bar.node2.position.z);

    freeNodeArray(&a);
    freeFrameBarArray(&fb);
    return 0;
}
