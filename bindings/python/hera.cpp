#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>


namespace py = pybind11;

void init_bt(py::module&);
void init_ws(py::module&);
void init_ws_geom(py::module&);


PYBIND11_MODULE(_hera, m)
{
    m.doc() = "Hera python bindings";
    init_bt(m);
    init_ws(m);
    init_ws_geom(m);
}
