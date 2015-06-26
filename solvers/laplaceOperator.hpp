/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef LAPLACE_OPERATOR_H
#define LAPLACE_OPERATOR_H

#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types2d;

void applyOperator(GlobalArrayD& first, const GlobalArrayD& other)
{
  assert(first.domain() == other.domain());

  assert((other.ghost_width()[0] >= 1 && other.ghost_width()[1] >= 1)
         || first.domain().group.nprocs == 1);

  other.iupdate();

  const Domain& dom = first.domain();
  const coords edge = other.ghost_width();
  coords_range inner = { dom.local().lower + edge, dom.local().upper - edge };
  double h = 1.0 / (double(dom.n[0]) + 1.0);

  AccessD b(first);
  const AccessD x(other);

  //inner region
  dom.for_each(inner,[h,&b,&x](coords ii)
  {
      coords left;  left[0] = ii[0]-1;  left[1] = ii[1];
      coords right; right[0] = ii[0]+1; right[1] = ii[1];
      coords below; below[0] = ii[0];   below[1] = ii[1]-1; below[2] = ii[2];
      coords above; above[0] = ii[0];   above[1] = ii[1]+1;
      b(ii) = (4.0 * x(ii) - x(left) - x(below) - x(right) - x(above)) / (h*h);
  });

  other.iupdate().wait();

  //outer region stencil
  auto stencil = [&dom,&b,&x,h](coords ii)
  {
    const coords left  = {{ii[0]-1, ii[1]}};
    const coords right = {{ii[0]+1, ii[1]}};
    const coords below = {{ii[0], ii[1]-1}};
    const coords above = {{ii[0], ii[1]+1}};

    b(ii) = 4.0 * x(ii);

    if(ii[0] > 0) b(ii) -= x(left);
    if(ii[1] > 0) b(ii) -= x(below);
    if(ii[0] < dom.n[0]-1) b(ii) -= x(right);
    if(ii[1] < dom.n[1]-1) b(ii) -= x(above);
    b(ii) /= h*h;

  };

  //treats the outer region (4 planes in 2D)

  for(int d = 0; d < 2; d++)
  {
      coords_range outer = { dom.local().lower, dom.local().upper };

      for(int di = 0; di < d; di++)
      {
          outer.lower[di] += edge[di];
          outer.upper[di] -= edge[di];
      }

      outer.upper[d] = dom.local().lower[d] + edge[d];
      dom.for_each(outer, stencil);
      outer.lower[d] = dom.local().upper[d] - edge[d];
      outer.upper[d] = dom.local().upper[d];
      dom.for_each(outer, stencil);
  }
}

#endif
