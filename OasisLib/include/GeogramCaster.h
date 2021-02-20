#ifndef GEOGRAMCASTER_H
#define GEOGRAMCASTER_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <geogram/basic/vecg.h>

namespace py = pybind11;

namespace pybind11
{
    namespace detail
    {
        template <index_t DIM, typename T>
        struct type_caster<vecng<DIM, T>>
        {
        public:
            #define COMMA ,
            PYBIND11_TYPE_CASTER(vecng<DIM COMMA T>, _("vecng<DIM, T>"));

            static py::handle cast(const vecng<DIM, T> &src, 
                    py::return_value_policy policy, py::handle parent)
            {
                py::array_t<T> a({DIM}, {sizeof(T)}, src.data() );
                return a.release();
            }
        };
    }
}

#endif