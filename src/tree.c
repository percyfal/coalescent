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
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
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
	printf(", mutations: %i, time: %.2f, visited: %i\n", n->mutations, n->time, n->visited);
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
 * reset_tree - mark all nodes as unvisited
 * 
 * @param n - Node
 */void reset_tree(struct Node *n)
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
	Rprintf("calculating segregating sites\n");
	while (n != NULL) {
		if (n != NULL) {
			Rprintf("id: %i, mutations: %i, visited: %i\n", n->id, n->mutations, n->visited);
			if (!n->visited) {
				Rprintf("Adding %i mutations from id %i\n", n->mutations, n->id);
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
	Rprintf("calculating tbl...\n");
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

struct Node *coalescent_tree(unsigned int N, double theta)
{
	// http://users.stat.umn.edu/~geyer/rc/
	GetRNGstate();
	int k;
	unsigned int i, j, n_nodes;
	// there are 2N-2 branches (rooted tree), which implies 2N-1 nodes
	n_nodes = 2*N - 1;
	struct Node *nodes[n_nodes];
	for (k=0; k < n_nodes; k++)
		nodes[k] = new_node(k);
	
	double rate, Ti, Ttot;
	Ttot = 0.0;
	i = j = N;
	while (i > 1) {
		// Set the coalescence time
		Rprintf("i: %i\n", i);
		rate = i * (i - 1) / 2.0;
		/* rate is in unit 1/s; want unit s so invert */
		Ti = rexp(1/rate);
		Ttot += Ti;
		nodes[j]->time = Ttot;

		// Make a number of mutations and sprinkle them out
		int n_mut, l;
		n_mut = (int) rpois((double) (theta / 2 * Ti));
		for (k=0; k<n_mut; k++) {
			l = (int) runif(0.0, (double) i);
			nodes[l]->mutations++;
		}
		Rprintf("Time (tau) %.2f, time %.2f, n mutations: %i\n", Ti, Ttot, n_mut);

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

	
SEXP tree(SEXP N, SEXP theta)
{
	SEXP ans;
	SEXP TMRCA = PROTECT(allocVector(REALSXP, 1));
	SEXP TBL = PROTECT(allocVector(REALSXP, 1));
	SEXP NMUT = PROTECT(allocVector(INTSXP, 1));
	struct Node *n;
	n = coalescent_tree(asInteger(N), asReal(theta));
	ans = PROTECT(allocVector(VECSXP, 3));
	NMUT = ScalarInteger(segregating_sites(n));
	TMRCA = ScalarReal(tmrca(n));
	TBL = ScalarReal(total_branch_length(n));
	SET_VECTOR_ELT(ans, 0, TMRCA);
	SET_VECTOR_ELT(ans, 1, TBL);
	SET_VECTOR_ELT(ans, 2, NMUT);
    UNPROTECT(4);
	return (ans);
}

/* int main(int argc, const char *argv) */
/* { */
/* 	time_t tt; */
/* 	tt = time(NULL); */
	
/* 	set_seed(tt,77911);   */
/* 	int k, n_nodes; */
/* 	unsigned int i, j, N; */
/* 	double theta; */
	
/* 	// samples */
/* 	N = 10; */
/* 	i = j = N; */
/* 	// Remember: E_TMRCA = 2(1-1/n), so in range [1,2). Theta = 1.0 */
/* 	// will then give 1-2 mutations on average which may be too low */
/* 	theta = 10.0; */
/* 	// there are 2N-2 branches (rooted tree), which implies 2N-1 nodes */
/* 	n_nodes = 2*N - 1; */
/* 	struct Node *nodes[n_nodes]; */
/* 	for (k=0; k < n_nodes; k++)  */
/* 		nodes[k] = new_node(k); */

/* 	double rate, Ti, Ttot; */
/* 	Ttot = 0.0; */
/* 	while (i > 1) { */
/* 		// Set the coalescence time */
/* 		rate = i * (i - 1) / 2.0; */
/* 		/\* rate is in unit 1/s; want unit s so invert *\/ */
/* 		Ti = rexp(1/rate); */
/* 		Ttot += Ti; */
/* 		nodes[j]->time = Ttot; */

/* 		// Make a number of mutations and sprinkle them out */
/* 		int n_mut, l; */
/* 		n_mut = (int) rpois((double) (theta / 2 * Ti)); */
/* 		for (k=0; k<n_mut; k++) { */
/* 			l = (int) runif(0.0, (double) i); */
/* 			nodes[l]->mutations++; */
/* 		} */
/* 		printf("Time (tau) %.2f, time %.2f, n mutations: %i\n", Ti, Ttot, n_mut); */

/* 		// Note that the coalescent event can be performed after */
/* 		// setting the time and placing mutations as we know the id of */
/* 		// the parent beforehand */
/* 		coalesce(nodes, j, i, theta); */
/* 		i--; */
/* 		j++; */
/* 	} */
/* 	for (k=0; k<n_nodes; k++) */
/* 		print_node(nodes[k]); */
/* 	for (k=0; k<n_nodes; k++) { */
/* 		if (isroot(nodes[k])) */
/* 			print_node(nodes[k]); */
/* 	} */
	

/* 	return 0; */
/* } */
