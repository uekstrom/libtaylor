// FIXME Proper include paths (needs libtaylor to be a proper CMake target)
#include <polymul.hpp>
#include <taylor.hpp>

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(pytaylor, m) {
  py::class_<polymul::polynomial<double, 1, 1>, std::shared_ptr<polymul::polynomial<double, 1, 1>>>(m, "polynomial")
    .def(py::init<double>());

  py::class_< taylor<double, 1, 1>, std::shared_ptr<taylor<double, 1, 1>>, polymul::polynomial<double, 1, 1>>(m, "taylor")
    .def(py::init<double>())
    .def("deriv_facs", &taylor<double, 1, 1>::deriv_facs);
}
