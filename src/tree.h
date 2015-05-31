/**
 * @file   tree.h
 * @author Per Unneberg
 * @date   Fri May 22 16:21:14 2015
 *
 * @brief   
 *
 *
 */

#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED

typedef int boolean;
#define TRUE 1
#define FALSE 0

// http://stackoverflow.com/questions/8257714/how-to-convert-an-int-to-string-in-c
// #define MAX_INT_CHARSIZE ((CHAR_BIT * sizeof(int) - 1) / 3 + 2)
#define NEWICK_BUFSIZE 10

struct Node {
	struct Node *parent;
	struct Node *left;
	struct Node *right;
	double time;
	unsigned int mutations;
	unsigned int id;
	boolean visited;
};


void print_node(struct Node *n);

#endif
