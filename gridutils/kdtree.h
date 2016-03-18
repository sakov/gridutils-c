/*
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2009 John Tsiombikas <nuclear@siggraph.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct kdtree kdtree;
typedef struct kdres kdres;
typedef struct kdnode kdnode;

/* create a kd-tree for "k"-dimensional data */
kdtree *kd_create(int k);

/* free the kdtree */
void kd_free(kdtree *tree);

/* remove all the elements from the tree */
void kd_clear(kdtree *tree);

/* if called with non-null 2nd argument, the function provided
 * will be called on data pointers (see kd_insert) when nodes
 * are to be removed from the tree.
 */
void kd_data_destructor(kdtree *tree, void (*destr)(void*));

/* insert a node, specifying its position, and optional data */
int kd_insert(kdtree *tree, const double *pos, void *data);

/* Find the nearest node from a given point.
 *
 * This function returns a pointer to a result set with at most one element.
 */
kdnode *kd_nearest(kdtree *tree, const double *pos);

/* Find any nearest nodes from a given point within a range.
 *
 * This function returns a pointer to a result set, which can be manipulated
 * by the kd_res_* functions.
 * The returned pointer can be null as an indication of an error. Otherwise
 * a valid result set is always returned which may contain 0 or more elements.
 * The result set must be deallocated with kd_res_free after use.
 */
kdres *kd_nearest_range(kdtree *tree, const double *pos, double range);

/* frees a result set returned by kd_nearest_range() */
void kd_res_free(kdres *set);

/* returns the size of the result set (in elements) */
int kd_res_size(kdres *set);

/* rewinds the result set iterator */
void kd_res_rewind(kdres *set);

/* returns non-zero if the set iterator reached the end after the last element */
int kd_res_end(kdres *set);

/* advances the result set iterator, returns non-zero on success, zero if
 * there are no more elements in the result set.
 */
int kd_res_next(kdres *set);

/* returns the data pointer (can be null) of the current result set item
 * and optionally sets its position to the pointers(s) if not null.
 */
void *kd_res_item(kdres *set, double *pos);

/* equivalent to kd_res_item(set, 0) */
void *kd_res_item_data(kdres *set);

void* kdnode_data(kdnode *node);
double* kdnode_pos(kdnode *node);

#ifdef __cplusplus
}
#endif

#endif	/* _KDTREE_H_ */
