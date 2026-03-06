#include <nano_geo_matrix/bessel/myBessel_bl.hpp>          // satisfies mie/geometry/single.hpp
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/mie/geometry/gimme_p.hpp>  // satisfies cup/mie/frohlich.hpp expectation of pcfc0, ..., pcfc3
#include <nano_geo_matrix/mie/geometry/single.hpp>      // satisfies CUP_HAS_FROHLICH_GEOMETRY expectation
#include <nano_geo_matrix/cup/cup.hpp>

int main() {
    nanosphere s;
    s.init();
    return 0;
}