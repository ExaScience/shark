#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types2d;

// Laplacian used in the CG, borrowed from Helsim
// To temporarily be used in place of the async version below (needs more work)
// Treats outer and inner regions differently:
//   the stencil in the inner region knows its neighbours exist
//   the stencil in the outer region needs to check

//BEGIN IMEN

void compute_diffusionConvection(GlobalArrayD& first, const GlobalArrayD& other)
{
	  assert(first.domain() == other.domain());

	  assert((other.ghost_width()[0] >= 1 && other.ghost_width()[1] >= 1)
	         || first.domain().group.nprocs == 1);

	  other.update();

	  const Domain& dom = first.domain();
	  const coords edge = {{1,1}};
	  coords_range inner = { dom.total().lower + edge, dom.total().upper - edge };
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
	      b(ii) = (x(ii) - x(left) - 3 * x(below) + 3 * x(right) + x(above)) / (h*h);
	  });

	  //outer region stencil
	  auto stencil = [&dom,&b,&x,h](coords ii)
	  {
		/*if j > 0: A[j*n+i, j*n+i-n] = -1
		  if i > 0: A[j*n+i, j*n+i-1] = -3
		  if i < n-1: A[j*n+i, j*n+i+1] = 3
		  if j < n-1: A[j*n+i, j*n+i+n] = 1
		 */

	    const coords left  = {{ii[0]-1, ii[1]}};
	    const coords right = {{ii[0]+1, ii[1]}};
	    const coords below = {{ii[0], ii[1]-1}};
	    const coords above = {{ii[0], ii[1]+1}};

	    b(ii) = x(ii);

	    if(ii[0] > 0) b(ii) -= x(left);
	    if(ii[1] > 0) b(ii) = b(ii) - 3 * x(below);
	    if(ii[0] < dom.n[0]-1) b(ii) = b(ii) + 3 * x(right);
	    if(ii[1] < dom.n[1]-1) b(ii) += x(above);

	    b(ii) /= h*h;

	  };

	  //treats the outer region (4 planes in 2D)

	  for(int d = 0; d < 2; d++)
	  {
	      coords_range outer = { dom.total().lower, dom.total().upper };

	      for(int di = 0; di < d; di++)
	      {
	          outer.lower[di] += edge[di];
	          outer.upper[di] -= edge[di];
	      }

	      outer.upper[d] = dom.total().lower[d] + edge[d];
	      dom.for_each(outer, stencil);
	      outer.lower[d] = dom.total().upper[d] - edge[d];
	      outer.upper[d] = dom.total().upper[d];
	  }
}

//END IMEN

void applyOperator(GlobalArrayD& first, const GlobalArrayD& other)
{
  assert(first.domain() == other.domain());

  assert((other.ghost_width()[0] >= 1 && other.ghost_width()[1] >= 1)
         || first.domain().group.nprocs == 1);

  other.update();

  const Domain& dom = first.domain();
  const coords edge = {{1,1}};
  coords_range inner = { dom.total().lower + edge, dom.total().upper - edge };
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
      coords_range outer = { dom.total().lower, dom.total().upper };

      for(int di = 0; di < d; di++)
      {
          outer.lower[di] += edge[di];
          outer.upper[di] -= edge[di];
      }

      outer.upper[d] = dom.total().lower[d] + edge[d];
      dom.for_each(outer, stencil);
      outer.lower[d] = dom.total().upper[d] - edge[d];
      outer.upper[d] = dom.total().upper[d];
      dom.for_each(outer, stencil);
  }
}

/*Roel: needs some more work - temporarily disabling and reverting to old, non-async version

void applyOperator(GlobalArrayD& first, const GlobalArrayD& other)
{
  assert(first.domain() == other.domain());

  // assert((other.ghost_width()[0] >= 1 && other.ghost_width()[1] >= 1 && other.ghost_width()[2] >= 1) 
  // || first.domain().group.nprocs == 1);

  assert((other.ghost_width()[0] >= 1 && other.ghost_width()[1] >= 1) || first.domain().group.nprocs == 1);

  const coords& n = first.domain().n;

  double h = 1.0 / (double(n[0]) + 1.0);

  auto stencil = [&n, h](AccessD& b, const AccessD& x, coords ii)
  {
    const coords left  = {{ii[0]-1, ii[1]}};
    const coords right = {{ii[0]+1, ii[1]}};
    const coords below = {{ii[0], ii[1]-1}};
    const coords above = {{ii[0], ii[1]+1}};
    
    b(ii) = 4.0 * x(ii);
    if(ii[0] > 0) b(ii) -= x(left);
    if(ii[1] > 0) b(ii) -= x(below);
    if(ii[0] < n[0]-1) b(ii) -= x(right);
    if(ii[1] < n[1]-1) b(ii) -= x(above);
    b(ii) /= h*h;
  };

  const coords gw = other.ghost_width();
  const coords_range local = first.domain().local();

  other.updateBegin();

  const coords_range inner = { local.lower + gw, local.upper - gw };

  Region(first.domain(), inner).map(first, other, stencil);

  other.updateWait();

  for(int d = 0; d < 2; d++)
  {
    coords_range outer = { local.lower, local.upper };

    for(int di = 0; di < d; di++)
    {
      outer.lower[di] += gw[di];
      outer.upper[di] -= gw[di];
    }

    outer.upper[d] = local.lower[d] + gw[d];

    Region(first.domain(), outer).map(first, other, stencil);

    outer.lower[d] = local.upper[d] - gw[d];

    outer.upper[d] = local.upper[d];

    Region(first.domain(), outer).map(first, other, stencil);
  }

}*/

