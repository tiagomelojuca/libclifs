#include <stdio.h>
#include "../src/libclifs.h"

int main()
{
    GlobalSystem g;
    initGlobalSystem(&g);

    Node n1, n2;
    n1.position.x = 5.0;
    n2.position.y = 7.0;
    
    insertNodeGlobalSystem(&g, n1);
    insertNodeGlobalSystem(&g, n2);

    printf("%f\n", g.nodeArray.nodes[1].position.y);
    printf("%i\n", g.nodeArray.used);
    printf("all fine");

    freeGlobalSystem(&g);
    return 0;
}
