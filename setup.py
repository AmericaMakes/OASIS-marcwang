import os
import sys
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cfg = "Debug" if self.debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DPYTHON_EXECUTABLE={}".format(sys.executable),
            "-DEXAMPLE_VERSION_INFO={}".format(self.distribution.get_version()),
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
        ]
        build_args = []

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += ["-j{}".format(self.parallel)]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        try:
            subprocess.check_call(
                ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
            )
            subprocess.check_call(
                ["cmake", "--build", "."] + build_args, cwd=self.build_temp
            )
        except subprocess.CalledProcessError as e:
            print(e.output)

# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="OasisLib",
    version="0.0.1",
    author="Marc Wang",
    author_email="marc.wang@uwaterloo.ca",
    description="A test project using pybind11 and CMake",
    long_description="",
    ext_modules=[CMakeExtension("OasisLib")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)