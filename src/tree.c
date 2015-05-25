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
/* #include <R.h> */
/* See http://stat-www.berkeley.edu/classes/s243/rmath.html for instructions on compiling standalone */

#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "tree.h"

typedef int boolean;
#define TRUE 1
#define FALSE 0

/* http://stackoverflow.com/questions/14768230/malloc-for-struct-and-pointer-in-c */
struct node *newNode (size_t id)
{
	struct node *retval = malloc (sizeof (struct node));
	if (retval == NULL)
		return NULL;
	retval->parent = malloc (sizeof (struct node));
	retval->left = malloc (sizeof (struct node));
	retval->right = malloc (sizeof (struct node));
	retval->parent = NULL;
	retval->left = NULL;
	retval->right = NULL;
	retval->mutations = 0;
	retval->time = 0.0;
	retval->id = id;
	return retval;
}
void delNode (struct node *n) 
{
	if (n != NULL) {
		free (n->parent);
		free (n->left);
		free (n->right);
		free (n);
	}
}

void printNode(struct node *n)
{
	printf("Node id: %zu", n->id);
	if (n->parent == NULL)
		printf(", parent id: %p", n->parent);
	else
		printf(", parent id: %zu", n->parent->id);
	if (n->left == NULL)
		printf(", left id: %p", n->left);
	else
		printf(", left id: %zu", n->left->id);
	if (n->right == NULL)
		printf(", right id: %p", n->right);
	else
		printf(", right id: %zu", n->right->id);
	printf(", mutations: %zu, time: %.2f\n", n->mutations, n->time);
	
}

/** 
 * isleaf - determine whether node is leaf or not
 * 
 * @param n 
 * 
 * @return boolean
 */
boolean isleaf(struct node *n) 
{
	return ((n->left == NULL) && (n->right == NULL)) ? TRUE : FALSE;
}
/** 
 * isroot - determine whether node is root or not
 * 
 * @param n 
 * 
 * @return boolean
 */
boolean isroot(struct node *n) 
{
	return (n->parent == NULL) ? TRUE : FALSE;
}


	
int main(int argc, const char *argv)
{
	time_t tt;
	tt = time(NULL);
	
/* the numbers passed to set_seed are arbitrary -- 
   use constants for a reproducible sequence       */
	
	set_seed(tt,77911);  
	struct node *n;
	int i,j;
	i = j = 3;
	size_t id;
	id = 1;
	double Ti, Ttot, rate;
	Ttot = 0.0;
	while (i > 1) {
		j++;
		rate = i * (i - 1) / 2;
		Ti = rexp(rate);
		Ttot += Ti;
		i--;
		printf("Now i: %.2f\n", Ti);
		
	}
	
}

		
