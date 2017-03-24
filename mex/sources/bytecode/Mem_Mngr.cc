/*
 * Copyright (C) 2007-2011 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Mem_Mngr.hh"

Mem_Mngr::Mem_Mngr()
{
  swp_f = false;
  swp_f_b = 0;
}
/*void
Mem_Mngr::Print_heap()
{
  unsigned int i;
  mexPrintf("i   :");
  for (i = 0; i < CHUNK_SIZE; i++)
    mexPrintf("%3d ", i);
  mexPrintf("\n");
}
*/

void
Mem_Mngr::init_Mem()
{
  Chunk_Stack.clear();
  CHUNK_SIZE = 0;
  Nb_CHUNK = 0;
  NZE_Mem = NULL;
  NZE_Mem_add = NULL;
  CHUNK_heap_pos = 0;
  NZE_Mem_Allocated.clear();
}

void
Mem_Mngr::fixe_file_name(string filename_arg)
{
  filename_mem = filename_arg;
}

void
Mem_Mngr::init_CHUNK_BLCK_SIZE(int u_count)
{
  CHUNK_BLCK_SIZE = u_count;
}

NonZeroElem *
Mem_Mngr::mxMalloc_NZE()
{
  long unsigned int i;
  if (!Chunk_Stack.empty())           /*An unused block of memory available inside the heap*/
    {
      NonZeroElem *p1 = Chunk_Stack.back();
      Chunk_Stack.pop_back();
      return (p1);
    }
  else if (CHUNK_heap_pos < CHUNK_SIZE) /*there is enough allocated memory space available we keep it at the top of the heap*/
    {
      i = CHUNK_heap_pos++;
      return (NZE_Mem_add[i]);
    }
  else                                /*We have to allocate extra memory space*/
    {
      CHUNK_SIZE += CHUNK_BLCK_SIZE;
      Nb_CHUNK++;
      NZE_Mem = (NonZeroElem *) mxMalloc(CHUNK_BLCK_SIZE*sizeof(NonZeroElem));      /*The block of memory allocated*/
      error_msg.test_mxMalloc(NZE_Mem, __LINE__, __FILE__, __func__, CHUNK_BLCK_SIZE*sizeof(NonZeroElem));
      NZE_Mem_Allocated.push_back(NZE_Mem);
      if (!NZE_Mem)
        mexPrintf("Not enough memory available\n");
      if (NZE_Mem_add)
        {
          NZE_Mem_add = (NonZeroElem **) mxRealloc(NZE_Mem_add, CHUNK_SIZE*sizeof(NonZeroElem *));                                                                                                     /*We have to redefine the size of pointer on the memory*/
          error_msg.test_mxMalloc(NZE_Mem_add , __LINE__, __FILE__, __func__, CHUNK_SIZE*sizeof(NonZeroElem *));
        }
      else
        {
          NZE_Mem_add = (NonZeroElem **) mxMalloc(CHUNK_SIZE*sizeof(NonZeroElem *));                                                                                       /*We have to define the size of pointer on the memory*/
          error_msg.test_mxMalloc(NZE_Mem_add , __LINE__, __FILE__, __func__, CHUNK_SIZE*sizeof(NonZeroElem *));
        }

      if (!NZE_Mem_add)
        mexPrintf("Not enough memory available\n");
      for (i = CHUNK_heap_pos; i < CHUNK_SIZE; i++)
        NZE_Mem_add[i] = (NonZeroElem *) (NZE_Mem+(i-CHUNK_heap_pos));
      i = CHUNK_heap_pos++;
      return (NZE_Mem_add[i]);
    }
}

void
Mem_Mngr::mxFree_NZE(void *pos)
{
  unsigned int i;
  size_t gap;
  for (i = 0; i < Nb_CHUNK; i++)
    {
      gap = ((size_t) (pos)-(size_t) (NZE_Mem_add[i*CHUNK_BLCK_SIZE]))/sizeof(NonZeroElem);
      if ((gap < CHUNK_BLCK_SIZE) && (gap >= 0))
        break;
    }
  Chunk_Stack.push_back((NonZeroElem *) pos);
}

void
Mem_Mngr::Free_All()
{
  while (NZE_Mem_Allocated.size())
    {
      mxFree(NZE_Mem_Allocated.back());
      NZE_Mem_Allocated.pop_back();
    }
  if (NZE_Mem_add)
    mxFree(NZE_Mem_add);
  init_Mem();
}
