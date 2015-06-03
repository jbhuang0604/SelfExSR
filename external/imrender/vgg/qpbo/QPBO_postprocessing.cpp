/* QPBO_postprocessing.cpp */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "QPBO.h"

template <typename REAL>
	void QPBO<REAL>::ComputeWeakPersistencies()
{
	if (stage == 0) return;

	Node* i;
	Node* j;
	Node* stack = NULL;
	int component;

	for (i=nodes[0]; i<node_last[0]; i++)
	{
		code_assert(i->label>=-1 && i->label<=1);

		Node* i1 = GetMate0(i);

		if (i->label >= 0)
		{
			i->dfs_parent = i;
			i1->dfs_parent = i1;
			i->region = i1->region = 0;
		}
		else
		{
			i->dfs_parent = i1->dfs_parent = NULL;
			i->region = i1->region = -1;
		}
	}

	// first DFS
	for (i=nodes[0]; i<node_last[1]; i++)
	{
		if (i == node_last[0]) i = nodes[1];
		if (i->dfs_parent) continue;

		// DFS starting from i
		i->dfs_parent = i;
		i->dfs_current = i->first;
		while ( 1 )
		{
			if (!i->dfs_current)
			{
				i->next = stack;
				stack = i;

				if (i->dfs_parent == i) break;
				i = i->dfs_parent;
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j = i->dfs_current->head;
			if (!(i->dfs_current->r_cap>0) || j->dfs_parent)
			{
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j->dfs_parent = i;
			i = j;
			i->dfs_current = i->first;
		}
	}

	// second DFS
	component = 0;
	while ( stack )
	{
		i = stack;
		stack = i->next;
		if (i->region > 0) continue;

		i->region = ++ component;
		i->dfs_parent = i;
		i->dfs_current = i->first;
		while ( 1 )
		{
			if (!i->dfs_current)
			{
				if (i->dfs_parent == i) break;
				i = i->dfs_parent;
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j = i->dfs_current->head;
			if (!(i->dfs_current->sister->r_cap>0) || j->region>=0)
			{
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j->dfs_parent = i;
			i = j;
			i->dfs_current = i->first;
			i->region = component;
		}
	}

	// assigning labels
	for (i=nodes[0]; i<node_last[0]; i++)
	{
		if (i->label < 0)
		{
			code_assert(i->region > 0);
			if      (i->region > GetMate0(i)->region) { i->label = 0; i->region = 0; }
			else if (i->region < GetMate0(i)->region) { i->label = 1; i->region = 0; }
		}
		else code_assert(i->region == 0);
	}
}

template <typename REAL>
	void QPBO<REAL>::Stitch()
{
	if (stage == 0) return;

	Node* i;
	Node* i_mate;
	Node* j;
	Arc* a;
	Arc* a_mate;

	for (a=arcs[0], a_mate=arcs[1]; a<arc_max[0]; a++, a_mate++)
	if (a->sister)
	{
		a->r_cap = a_mate->r_cap = a->r_cap + a_mate->r_cap;

		i = a->sister->head;
		j = a->head;

		if (i->region==0 || i->region != j->region) continue;
		if (IsNode0(i))
		{
			if (i->user_label != 0) continue;
		}
		else
		{
			if (GetMate1(i)->user_label != 1) continue;
		}
		if (IsNode0(j))
		{
			if (j->user_label != 1) continue;
		}
		else
		{
			if (GetMate1(j)->user_label != 0) continue;
		}

		a->r_cap = a_mate->r_cap = 0;
	}

	for (i=nodes[0], i_mate=nodes[1]; i<node_last[0]; i++, i_mate++)
	{
		i->tr_cap = i->tr_cap - i_mate->tr_cap;
		i_mate->tr_cap = -i->tr_cap;
	}

	ComputeWeakPersistencies();
}

#include "instances.inc"
