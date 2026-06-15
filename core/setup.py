from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

ext_modules = [
    Pybind11Extension(
        "omp_core",
        sources=[
            "src/bindings.cpp",
            "src/tle.cpp",
            "src/sgp4.cpp",
            "src/twobody.cpp",
            "src/transforms.cpp",
            "src/transfer.cpp",
        ],
        include_dirs=["include"],
        cxx_std=17,
    ),
]

setup(
    name="omp_core",
    version="0.1.0",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
