#include <boost/python.hpp>
#include "dejavu_api.h"

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

BOOST_PYTHON_MODULE(libdejavu_python) {
        using namespace boost::python;
        def("random_paths", random_paths);
}
