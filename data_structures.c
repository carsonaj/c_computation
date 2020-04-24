
// linked list node

void init_Lnode(Lnode *node) {
    node->data = 0;
    node->next = NULL;
}
//----------------------------------------------------------------------------

//binary tree node

void init_Bnode(Bnode *node) {
    node->data = 0;
    node->left = NULL;
    node->right = NULL;
}

void append_Bnode(Bnode *child, Bnode *parent, int side) {
    if (side == -1) {
        if (parent->left == NULL) {
            parent->left = child;
            child->parent = parent;
        }
    }
    else if (side == 1) {
        if (parent->right == NULL) {
            parent->right = child;
            child->parent = parent;
        }
    }

}
