#include <stdlib.h>

// linked list node

typedef struct _Lnode Lnode;

struct _Lnode {
    double data;
    struct _Lnode *next;
};

void init_Lnode(Lnode *node);
//----------------------------------------------------------------------------

//binary tree node

typedef struct _Bnode Bnode;

struct _Bnode {
    int data;
    struct _Bnode *parent;
    struct _Bnode *left;
    struct _Bnode *right;
};

void init_Bnode(Bnode *node);
void append_Bnode(Bnode *child, Bnode *parent, int side);
//----------------------------------------------------------------------------
