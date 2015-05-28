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

struct Node {
	struct Node *parent;
	struct Node *left;
	struct Node *right;
	double time;
	unsigned int mutations;
	unsigned int id;
	boolean visited;
	char *state;
};


void print_node(struct Node *n);

#endif
