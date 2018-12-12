#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <memory>

// FIXME Proper include paths (needs libtaylor to be a proper CMake target)
#include <polymul.hpp>
#include <taylor.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

namespace {
template <typename T, int Nvar, int Ndeg>
void declarePolynomial(py::module & mod, const std::string & suffix) {
  using Class = polymul::polynomial<T, Nvar, Ndeg>;
  using PyClass = py::class_<Class, std::shared_ptr<Class>>;

  PyClass cls(mod, ("polynomial" + suffix).c_str());

  cls.def(py::init<T>());
}

template <typename T, int Nvar, int Ndeg>
void declareTaylor(py::module & mod, const std::string & suffix) {
  using Class = taylor<T, Nvar, Ndeg>;
  using PyClass =
      py::class_<Class, std::shared_ptr<Class>, polymul::polynomial<T, Nvar, Ndeg>>;

  PyClass cls(mod, ("taylor" + suffix).c_str());

  cls.def(py::init<T>());
  cls.def(py::init<T, int>());
  cls.def(py::init<T, int, T>());
  cls.def("deriv_facs", &Class::deriv_facs);
  cls.def(
      "__getitem__", py::overload_cast<int>(&Class::operator[]), py::is_operator());
  cls.def("__getitem__",
          py::overload_cast<int>(&Class::operator[], py::const_),
          py::is_operator());
  // TODO operators need to be revisited
  cls.def("__add__",
          [](const Class & self, const Class & other) { return self + other; },
          py::is_operator());
  cls.def("__sub__",
          [](const Class & self, const Class & other) { return self - other; },
          py::is_operator());
  cls.def("__mul__",
          [](const Class & self, const Class & other) { return self * other; },
          py::is_operator());
  cls.def("__truediv__",
          [](const Class & self, const Class & other) { return self / other; },
          py::is_operator());
}
} // namespace

PYBIND11_MODULE(_pytaylor, mod) {
  mod.doc() = "PyTaylor: Forward-mode Automatic Differentiation";

  declarePolynomial<double, 1, 1>(mod, "D_1_1");
  declareTaylor<double, 1, 1>(mod, "D_1_1");

  declarePolynomial<double, 2, 1>(mod, "D_2_1");
  declareTaylor<double, 2, 1>(mod, "D_2_1");
}
