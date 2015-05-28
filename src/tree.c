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
// #define MATHLIB_STANDALONE
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"

#include "tree.h"

/* http://stackoverflow.com/questions/14768230/malloc-for-struct-and-pointer-in-c */
/** 
 * new_node - create a new Node
 * 
 * @param id 
 * 
 * @return Node
 */
struct Node *new_node (unsigned int id, unsigned int L)
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
	// here we go with chars again - there must be an easier way to do
	// this
	retval->state = malloc(sizeof(char)*(L+1));
	char state [L];
	strcpy(state, "");
	for (int i=0; i<L; i++)
		strcat(state, "0");
	strcpy(retval->state, state);
	return retval;
}
/** 
 * delete_node - delete a Node
 * 
 * @param n 
 */
void delete_node (struct Node *n) 
{
	if (n != NULL) {
		// freeing n->parent from a leaf apparently 
		// frees the parent as n->parent itself is pointing to a
		// pointer? ... so don't do this. 
		//free (n->parent);
		free (n->left);
		free (n->right);
		free (n->state);
		free (n);
	}
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
/** 
 * remove_tree - free up memory
 *
 * @param n - Node
 */
void remove_tree(struct Node *n)
{
	int ndel = 0;
	boolean left = FALSE;
	while (n != NULL) {
		if (isleaf(n)) {
			if (isroot(n)) {
				ndel += 1;
				delete_node(n);
				n = NULL;
			} else {
				left = TRUE;
				if (n->parent->left == NULL)
					left = FALSE;
				else if (n->parent->left->id != n->id)
					left = FALSE;
				n = n->parent;
				(left == TRUE) ? delete_node(n->left) : delete_node(n->right);
				if (left == TRUE) 
					n->left = NULL;
				else
					n->right = NULL;
				ndel += 1;
			}
		} else if (n->left != NULL) {
			n = n->left;
		} else if (n->right != NULL) {
			n = n->right;
		} else {
			// Add error message?
			n = NULL;
		}
	}
	// debug macro wanted
	// fprintf(stderr, "deleted %i nodes\n", ndel);
}

/** 
 * visit_node - a depth-first traversal of the tree
 * 
 * @param n 
 * 
 * @return 
 */
struct Node *visit_node (struct Node *n)
{
	n->visited = TRUE;
	if (n->left != NULL && !n->left->visited)
		return n->left;
	else if (n->right != NULL && !n->right->visited)
		return n->right;
	else if (n->parent != NULL) {
		return n->parent;
	}
	else {
		//Rprintf("all nodes have been visited; returning null\n");
		return NULL;
	}
}
/** 
 * unvisit_node - mark a node as unvisited
 * 
 * @param n 
 * 
 * @return 
 */
struct Node *unvisit_node (struct Node *n)
{
	n->visited = FALSE;
	if (n->left != NULL && n->left->visited)
		return n->left;
	else if (n->right != NULL && n->right->visited)
		return n->right;
	else if (n->parent != NULL) {
		return n->parent;
	}
	else {
		//Rprintf("all nodes have been unvisited; returning null\n");
		return NULL;
	}
}
/** 
 * reset_tree - mark all nodes as unvisited
 * 
 * @param n - Node
 */
void reset_tree(struct Node *n)
{
	while (n != NULL) {
		if (n != NULL) {
			n->visited = FALSE;
			n = unvisit_node(n);
		}
	}
}
/** 
 * segregating_sites - calculate the number of segregating sites
 * 
 * @param n - Node
 * 
 * @return int - the number of segregating sites
 */
int segregating_sites(struct Node *n) 
{
	int sites = 0;
	reset_tree(n);
	//Rprintf("calculating segregating sites\n");
	while (n != NULL) {
		if (n != NULL) {
			//Rprintf("id: %i, mutations: %i, visited: %i\n", n->id, n->mutations, n->visited);
			if (!n->visited) {
				//Rprintf("Adding %i mutations from id %i\n", n->mutations, n->id);
				n->visited = TRUE;
				sites += n->mutations;
			}
		}
		n = visit_node(n);
	}
	return sites;
}
/** 
 * total_branch_length - calculate the total branch length. As each
 * node has recorded a coalescent event at the *total* time point, we
 * need to add the *differences* in time between node and its parent
 * 
 * @param n - Node
 * 
 * @return double - total branch length
 */
double total_branch_length(struct Node *n)
{
	double tbl = 0.0;
	reset_tree(n);
	//Rprintf("calculating tbl...\n");
	while (n != NULL) {
		if (n != NULL) {
			//Rprintf("id: %i, tbl: %.3f, visited: %i\n", n->id, n->time, n->visited);
			if (!n->visited && n->parent != NULL) {
				//Rprintf("Adding time %.3f  from id %i\n", n->parent->time - n->time, n->id);
				n->visited = TRUE;
				tbl += n->parent->time - n->time;
			}
		}
		n = visit_node(n);
	}
	return tbl;
}
/** 
 * tmrca - calculate time to most recent common ancestor. This is
 * simply the time recorded in the parent node
 * 
 * @param n 
 * 
 * @return 
 */
double tmrca(struct Node *n)
{
	if (isroot(n))
		return n->time;
	else
		printf ("ERROR: %i is not root!\n", n->id);
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
static void coalesce (struct Node* nodes[], unsigned int node_census, unsigned int current_sample_size, double theta)
{
	struct Node *tmp;
	int i, j, parent_id, left, right;
	parent_id = node_census;

	// Select left child
	left = (int) runif(0.0, (double) current_sample_size);
	// swap left with last element in current_sample
	tmp = nodes[current_sample_size - 1];
	nodes[current_sample_size - 1] = nodes[left];
	nodes[left] = tmp;
	
	// Select right child, sampling from current_sample_size - 1 nodes
	right = (int) runif(0.0, (double) (current_sample_size - 1));
	// swap right with parent
	tmp = nodes[parent_id];
	nodes[parent_id] = nodes[right];
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
}

/** 
 * 
 * 
 * @param N number of samples
 * @param theta mutation rate
 * @param L number of sites
 * 
 * @return root node
 */
struct Node *coalescent_tree(unsigned int N, double theta, unsigned int L)
{
	// http://users.stat.umn.edu/~geyer/rc/
	GetRNGstate();
	int k;
	unsigned int i, j, n_nodes;
	// there are 2N-2 branches (rooted tree), which implies 2N-1 nodes
	n_nodes = 2*N - 1;
	struct Node *nodes[n_nodes];
	for (k=0; k < n_nodes; k++)
		nodes[k] = new_node(k, L);
	
	double rate, Ti, Ttot;
	Ttot = 0.0;
	i = j = N;
	while (i > 1) {
		// Set the coalescence time
		rate = i * (i - 1) / 2.0;
		/* rate is in unit 1/s; want unit s so invert */
		Ti = rexp(1/rate);
		Ttot += Ti;
		nodes[j]->time = Ttot;

		// Make a number of mutations and sprinkle them out
		int n_mut, l;
		// Remember; we need to multiply rate also by the number of
		// branches; otherwise we're only generating a number of
		// mutations proportional to TMRCA, not TBL
		n_mut = (int) rpois((double) (theta / 2 * Ti * i));
		for (k=0; k<n_mut; k++) {
			l = (int) runif(0.0, (double) i);
			nodes[l]->mutations++;
		}

		// Note that the coalescent event can be performed after
		// setting the time and placing mutations as we know the id of
		// the parent beforehand
		coalesce(nodes, j, i, theta);
		i--;
		j++;
	}
	PutRNGstate();
	for (k=0; k < n_nodes; k++) {
		if (isroot(nodes[k]))
			return nodes[k];
	}
}

/** 
 * tree - get the coalescent
 * 
 * @param N - sample size 
 * @param theta - mutation rate
 * @param L - sequence length
 * 
 * @return 
 */	
SEXP tree(SEXP N, SEXP theta, SEXP L)
{
	SEXP ans;
	SEXP TMRCA = PROTECT(allocVector(REALSXP, 1));
	SEXP TBL = PROTECT(allocVector(REALSXP, 1));
	SEXP NMUT = PROTECT(allocVector(INTSXP, 1));
	struct Node *n;
	n = coalescent_tree(asInteger(N), asReal(theta), asInteger(L));
	ans = PROTECT(allocVector(VECSXP, 3));
	NMUT = ScalarInteger(segregating_sites(n));
	TMRCA = ScalarReal(tmrca(n));
	TBL = ScalarReal(total_branch_length(n));
	SET_VECTOR_ELT(ans, 0, TMRCA);
	SET_VECTOR_ELT(ans, 1, TBL);
	SET_VECTOR_ELT(ans, 2, NMUT);
    UNPROTECT(4);
	// Clean up memory
	remove_tree(n);
	return (ans);
}
/** 
 * print_node - print representation of Node
 * 
 * @param n 
 */
void print_node(struct Node *n)
{
	printf("Node id: %i", n->id);
	(n->parent == NULL) ? printf(", parent id: %p", n->parent) : printf(", parent id: %i", n->parent->id);
	(n->left == NULL) ? printf(", left id: %p", n->left) : printf(", left id: %i", n->left->id);
	(n->right == NULL) ? printf(", right id: %p", n->right) : printf(", right id: %i", n->right->id);
	printf(", mutations: %i, time: %.2f, visited: %i\t", n->mutations, n->time, n->visited);
	printf("%s\n", n->state);
}
