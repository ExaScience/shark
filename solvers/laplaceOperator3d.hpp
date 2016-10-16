#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types3d;


void applyOperator(GlobalArrayD& first, const GlobalArrayD& other, int degree)
{
    assert(first.domain() == other.domain());
        assert((other.ghost_width()[0] >= 1 && other.ghost_width()[1] >= 1 && other.ghost_width()[2] >= 1)
                || first.domain().group.nprocs == 1);

    const Domain& dom = other.domain();
    const coords edge = other.ghost_width();
    coords_range inner = { dom.local().lower + edge, dom.local().upper - edge };

    AccessD b(first);
    const AccessD x(other);

    auto f = other.iupdate();

    for (int r=0; r<degree; ++r) {
        f.test();
   
        //inner region
        dom.for_each(inner,[&b,&x](coords ii)
                {
                coords left;   left[0] = ii[0]-1;  left[1] = ii[1];   left[2] = ii[2];
                coords right; right[0] = ii[0]+1; right[1] = ii[1];   right[2] = ii[2];
                coords below; below[0] = ii[0];   below[1] = ii[1]-1; below[2] = ii[2];
                coords above; above[0] = ii[0];   above[1] = ii[1]+1; above[2] = ii[2];
                coords back;   back[0] = ii[0];    back[1] = ii[1];    back[2] = ii[2]-1;
                coords front; front[0] = ii[0];   front[1] = ii[1];   front[2] = ii[2]+1;

                b(ii) = 6.0 * x(ii) - x(left) - x(below) - x(above) - x(right) - x(back) - x(front);
                });
    }

    //outer region stencil
    f.wait();

    auto stencil = [&dom,&b,&x](coords ii)
    {
        const coords left  = {{ii[0]-1, ii[1], ii[2]}};
        const coords right = {{ii[0]+1, ii[1] , ii[2]}};
        const coords below = {{ii[0], ii[1]-1, ii[2] }};
        const coords above = {{ii[0], ii[1]+1, ii[2]}};
        const coords front = {{ii[0], ii[1], ii[2]-1 }};
        const coords back = {{ii[0], ii[1], ii[2]+1}};

        b(ii) = 6.0 * x(ii);

        if(ii[0] > 0)  b(ii) -= x(left);
        if(ii[1] > 0)  b(ii) -= x(below);
        if(ii[2] > 0)  b(ii) -= x(front);

        if(ii[0] < dom.n[0]-1)  b(ii) -= x(right);
        if(ii[1] < dom.n[1]-1)  b(ii) -= x(above);
        if(ii[2] < dom.n[2]-1)  b(ii) -= x(back);

    };

    //treats the outer region
    for (int r=0; r<degree; ++r) {
        for(int d = 0; d < 3; d++)
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
}

