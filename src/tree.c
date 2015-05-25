/**
 * @file   tree.c
 * @author Per Unneberg
 * @date   Fri May 22 16:18:23 2015
 *
 * @brief   
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
/* See http://stat-www.berkeley.edu/classes/s243/rmath.html for instructions on compiling standalone */
#define MATHLIB_STANDALONE
#include "Rmath.h"
#include "tree.h"

/* http://stackoverflow.com/questions/14768230/malloc-for-struct-and-pointer-in-c */
/** 
 * new_node - create a new Node
 * 
 * @param id 
 * 
 * @return Node
 */
struct Node *new_node (unsigned int id)
{
	struct Node *retval = malloc (sizeof (struct Node));
	if (retval == NULL)
		return NULL;
	retval->parent = malloc (sizeof (struct Node));
	retval->left = malloc (sizeof (struct Node));
	retval->right = malloc (sizeof (struct Node));
	retval->parent = NULL;
	retval->left = NULL;
	retval->right = NULL;
	retval->mutations = 0;
	retval->time = 0.0;
	retval->id = id;
	retval->visited = FALSE;
	return retval;
}
/** 
 * del_node - delete a Node
 * 
 * @param n 
 */
void del_node (struct Node *n) 
{
	if (n != NULL) {
		free (n->parent);
		free (n->left);
		free (n->right);
		free (n);
	}
}

boolean visited(struct Node *n) 
{
	return n->visited;
}
struct Node *descend_node (struct Node *n)
{
	if (n->left != NULL && !n->left->visited)
		return n->left;
	else if (n->right != NULL && !n->right->visited)
		return n->right;
	else {
		n->visited = TRUE;
		return n->parent;
	}
}

/** 
 * print_node - print representation of Node
 * 
 * @param n 
 */
void print_node(struct Node *n)
{
	printf("Node id: %i", n->id);
	if (n->parent == NULL)
		printf(", parent id: %p", n->parent);
	else
		printf(", parent id: %i", n->parent->id);
	if (n->left == NULL)
		printf(", left id: %p", n->left);
	else
		printf(", left id: %i", n->left->id);
	if (n->right == NULL)
		printf(", right id: %p", n->right);
	else
		printf(", right id: %i", n->right->id);
	printf(", mutations: %i, time: %.2f\n", n->mutations, n->time);
}
/** 
 * coalesce - perform a coalescent event
 *
 * Two node indices are sampled from the current_sample_size. The node
 * at index node_census + 1
 * 
 * @param nodes - an array of all nodes in tree
 *
 * @param node_census - the highest node index that has been assigned
 *                      children/parents (initialized as number of leaves)
 * @param current_sample_size - the number of nodes to currently sample from
 * 
 */
void coalesce (struct Node* nodes[], unsigned int node_census, unsigned int current_sample_size)
{
	struct Node *tmp;
	int parent_id, left, right;

	// Select left child
	left = (int) runif(0.0, (double) current_sample_size);
	// swap left with last element in current_sample
	tmp = nodes[current_sample_size - 1];
	nodes[current_sample_size - 1] = nodes[left];
	nodes[left] = tmp;
	
	// Select right child
	right = (int) runif(0.0, (double) (current_sample_size - 1));
	// swap right with last element in node_census + 1 which is the new parent
	tmp = nodes[node_census];
	nodes[node_census] = nodes[right];
	nodes[right] = tmp;

	// reset indices
	parent_id = right;
	left = current_sample_size - 1;
	right = node_census;
	
	// Set parent to left and child to parent
	nodes[left]->parent = nodes[parent_id];
	nodes[parent_id]->left = nodes[left];
	// Set parent to right and child to parent
	nodes[right]->parent = nodes[parent_id];
	nodes[parent_id]->right = nodes[right];

	// Set the coalescence time
	int i;
	double rate, Ti;
	i = current_sample_size;
	rate = i * (i - 1) / 2.0;
	/* rate is in unit 1/s; want unit s so invert */
	Ti = rexp(1/rate);
	//printf("Rate: %.2f, time: %.2f\n", rate, Ti);
	nodes[parent_id]->time = Ti;
}

/** 
 * isleaf - determine whether Node is leaf or not
 * 
 * @param n 
 * 
 * @return boolean
 */
boolean isleaf(struct Node *n) 
{
	return ((n->left == NULL) && (n->right == NULL)) ? TRUE : FALSE;
}
/** 
 * isroot - determine whether Node is root or not
 * 
 * @param n 
 * 
 * @return boolean
 */
boolean isroot(struct Node *n) 
{
	return (n->parent == NULL) ? TRUE : FALSE;
}

int main(int argc, const char *argv)
{
	time_t tt;
	tt = time(NULL);
	
	set_seed(tt,77911);  
	int k;
	unsigned int i, j, N;
	N = 100;
	i = j = N;
	// there are 2N-2 branches, so 2N-1 nodes
	int n_nodes;
	n_nodes = 2*N - 1;
	struct Node *nodes[n_nodes];
	for (k=0; k < n_nodes; k++) 
		nodes[k] = new_node(k);
	
	while (i > 1) {
		coalesce(nodes, j, i);
		i--;
		j++;
	}
	for (k=0; k<n_nodes; k++)
		print_node(nodes[k]);

	return 0;
}



		
