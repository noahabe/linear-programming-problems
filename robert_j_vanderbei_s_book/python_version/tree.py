#!/usr/bin/env python
from tnode import tnode
""" generated source for module tree """
# /
#     Copyright (c) Robert J. Vanderbei, 1994
#     All Rights Reserved
# /

""" generated source for class tree """
# public static final int TRUE = 1;
# public static final int FALSE = 0;

root = None
curnode = None

def killtree():
    """ generated source for method killtree """
    global root
    global curnode
    if root != None:
        killnode(root)
    root = None

def addtree(data):
    global root
    global curnode
    """ generated source for method addtree """
    node = root
    parent = None
    if root == None:
        root = tnode()
        # if (((root) = (tnode)malloc((1) * sizeof(tnode))) == null && (1) > 0)
        root.data = data
        root.parent = None
        root.left = None
        root.right = None
        return
    while node != None:
        if node.data > data:
            parent = node
            node = node.left
        elif node.data < data:
            parent = node
            node = node.right
        else:
            return
        #  already there
    node = tnode()
    # if (((node) = (tnode)malloc((1) * sizeof(tnode))) == null && (1) > 0)
    node.data = data
    node.parent = parent
    node.left = None
    node.right = None
    if parent.data > data:
        parent.left = node
    else:
        parent.right = node

def deltree(data):
    global root
    global curnode
    """ generated source for method deltree """
    node = root
    parent = None
    node1 = tnode()
    node2 = tnode()
    while node != None and node.data != data:
        if node.data > data:
            node = node.left
        elif node.data < data:
            node = node.right
    if node == None:
        return
    #  not there anyway
    parent = node.parent
    if node.right == None:
        node1 = node.left
        if node1 != None:
            node1.parent = parent
    elif node.left == None:
        node1 = node.right
        node1.parent = parent
    else:
        node1 = node.left
        node2 = node1
        while node2.right != None:
            node2 = node2.right
        node2.right = node.right
        node.right.parent = node2
        node.left.parent = parent
    if parent == None:
        root = node1
    else:
        if parent.data > data:
            parent.left = node1
        else:
            parent.right = node1
        #if (node) != None:
            # free((node));
        #(node) = None

def getfirst():
    global root
    global curnode
    """ generated source for method getfirst """
    node = None
    parent = None
    if root == None:
        return -1
    node = root
    while node != None:
        parent = node
        node = node.left
    curnode = parent
    return curnode.data



def getnext():
    global root
    global curnode
    """ generated source for method getnext """
    node = None
    par = None
    if curnode == None:
        return -1
    if curnode.right != None:
        node = curnode.right
        while node != None:
            par = node
            node = node.left
        curnode = par
        return curnode.data
    node = curnode.parent
    while node != None:
        if node.data > curnode.data:
            break
        node = node.parent
    curnode = node
    if curnode != None:
        return curnode.data
    else:
        return -1
        #  no more

def getlast():
    global root
    global curnode
    """ generated source for method getlast """
    node = None
    parent = None
    if root == None:
        return -1
    node = root
    while node != None:
        parent = node
        node = node.right
    curnode = parent
    return curnode.data

def getprev():
    global root
    global curnode
    """ generated source for method getprev """
    node = None
    par = None
    if curnode == None:
        return -1
    if curnode.left != None:
        node = curnode.left
        while node != None:
            par = node
            node = node.right
        curnode = par
        return curnode.data
    node = curnode.parent
    while node != None:
        if node.data < curnode.data:
            break
        node = node.parent
    curnode = node
    if curnode != None:
        return curnode.data
    else:
        return -1
        #  no more

def printtree():
    global root
    global curnode
    """ generated source for method printtree """
    if root != None:
        print('root node is {:4d}'.format(root.data))
    print("  data parent left right")
    printnode(root)


def killnode( node):
    """ generated source for method killnode """
    if node.left != None:
        killnode(node.left)
    if node.right != None:
        killnode(node.right)
    if node != None:
        #    if (node) != None:
        node = None


def printnode(node):
    """ generated source for method printnode """
    parent = int()
    left = int()
    right = int()
    if node != None:
        if node.parent != None:
            parent = node.parent.data
        else:
            parent = -1
        if node.left != None:
            left = node.left.data
        else:
            left = -1
        if node.right != None:
            right = node.right.data
        else:
            right = -1
        printnode(node.left)
        print("    {:4d} {:4d} {:4d} {:4d}".format(node.data, parent, left, right))
        printnode(node.right)

    #
    # 		   TEST PROGRAM
    # 	#define N 50
    # 	main()
    # 	{
    # 	    int data, outdata;
    # 	    int list[N], list2[N];
    # 	    int i,j,k;
    # 	    for (k=0; k<5; k++) {
    # 	        for (i=0; i<N; i++) {
    # 		    list[i]  = 0;
    # 		    list2[i] = 0;
    # 	        }
    # 	        for (i=0; i<3; i++) {
    # 		    data = N*drand48();
    # 		    printf("          %4d \n", data);
    # 		    list[data] = 1;
    # 		    addtree( data );
    # 	        }
    # 	        for (outdata = getfirst(); ; outdata=getnext()) {
    # 		    list2[outdata] = 1;
    # 	            printf("%4d \n", outdata);
    # 		    if (outdata >= N-1) break;
    # 	            for (i=0; i<3; i++) {
    # 		        data = outdata + 1 + (N-outdata-1)*drand48();
    # 		        printf("          %4d \n", data);
    # 		        list[data] = 1;
    # 		        addtree( data );
    # 	            }
    # 	        }
    # 		printtree();
    # 	        killtree();
    # 	        for (i=0; i<N; i++) {
    # 		    printf(" %4d %4d %4d ", i, list[i], list2[i]);
    # 		    if (list[i] == list2[i]) { printf("\n"); }
    # 		    else                     { printf(" * \n"); }
    # 	        }
    # 	    }
    # 	}
    #

