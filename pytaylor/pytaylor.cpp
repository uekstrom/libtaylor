#include <pybind11/pybind11.h>

// FIXME Proper include paths (needs libtaylor to be a proper CMake target)
#include <polymul.hpp>
#include <taylor.hpp>

namespace py = pybind11;

namespace {
template <typename T, int Nvar, int Ndeg>
void declarePolynomial(py::module & mod, const std::string & suffix) {
  using Class = polymul::polynomial<T, Nvar, Ndeg>;
  py::class_<Class, std::shared_ptr<Class>> cls(mod, ("polynomial" + suffix).c_str());

  cls.def(py::init<T>());
}

template <typename T, int Nvar, int Ndeg>
void declareTaylor(py::module & mod, const std::string & suffix) {
  using Class = taylor<T, Nvar, Ndeg>;
  py::class_<Class, std::shared_ptr<Class>, polymul::polynomial<T, Nvar, Ndeg>> cls(
      mod, ("taylor" + suffix).c_str());

  cls.def(py::init<T>());
  cls.def("deriv_facs", &Class::deriv_facs);
}
} // namespace

PYBIND11_MODULE(pytaylor, mod) {
  declarePolynomial<double, 1, 1>(mod, "D_1_1");
  declareTaylor<double, 1, 1>(mod, "D_1_1");
}
