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


struct node {
	struct node *parent;
	struct node *left;
	struct node *right;
	double time;
	size_t mutations;
	size_t id;
};

#endif
