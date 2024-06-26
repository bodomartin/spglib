import os
import sys
import sysconfig

setup_type = sys.argv[1]

try:
    from setuptools import Extension, setup
    from setuptools.command.build_ext import build_ext

    use_setuptools = True
    print("setuptools is used.")

    # Trick to pip install numpy when it's not installed.
    # Reference: https://stackoverflow.com/questions/2379898/
    class CustomBuildExtCommand(build_ext):
        def run(self):
            import numpy

            self.include_dirs.append(numpy.get_include())
            build_ext.run(self)

except ImportError:
    from distutils.core import Extension, setup

    use_setuptools = False
    print("distutils is used.")

    try:
        from numpy.distutils.misc_util import get_numpy_include_dirs
    except ImportError:
        print(
            "numpy.distutils.misc_util cannot be imported. Please install "
            "numpy first before installing spglib..."
        )
        sys.exit(1)

# Workaround Python issue 21121
config_var = sysconfig.get_config_var("CFLAGS")
if config_var is not None and "-Werror=declaration-after-statement" in config_var:
    os.environ["CFLAGS"] = config_var.replace("-Werror=declaration-after-statement", "")

sources = [
    "arithmetic.c",
    "cell.c",
    "delaunay.c",
    "debug.c",
    "determination.c",
    "hall_symbol.c",
    "kgrid.c",
    "kpoint.c",
    "magnetic_spacegroup.c",
    "mathfunc.c",
    "msg_database.c",
    "niggli.c",
    "overlap.c",
    "pointgroup.c",
    "primitive.c",
    "refinement.c",
    "sitesym_database.c",
    "site_symmetry.c",
    "spacegroup.c",
    "spin.c",
    "spg_database.c",
    "spglib.c",
    "symmetry.c",
]

if os.path.exists("src"):
    source_dir = "src"
else:
    source_dir = os.path.join(os.pardir, "src")

include_dirs = [
    source_dir,
]
if not use_setuptools:
    include_dirs += get_numpy_include_dirs()

for i, s in enumerate(sources):
    sources[i] = os.path.join(source_dir, s)

extra_compile_args = []
if setup_type == "test":
    extra_compile_args.append("-UNDEBUG")
extra_link_args = []
define_macros = []

# Uncomment to activate OpenMP support for gcc
# extra_compile_args += ['-fopenmp']
# extra_link_args += ['-lgomp']

# For debugging
# define_macros = [('SPGWARNING', None),
#                  ('SPGDEBUG', None)]

extension = Extension(
    "spglib._spglib",
    include_dirs=include_dirs,
    sources=["_spglib.c"] + sources,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    define_macros=define_macros,
)

version_nums = [None, None, None]
with open(os.path.join(source_dir, "version.h")) as w:
    for line in w:
        for i, chars in enumerate(("MAJOR", "MINOR", "MICRO")):
            if chars in line:
                version_nums[i] = int(line.split()[2])

# To deploy to pypi by travis-CI
nanoversion = 0
if os.path.isfile("__nanoversion__.txt"):
    with open("__nanoversion__.txt") as nv:
        try:
            for line in nv:
                nanoversion = int(line.strip())
                break
        except ValueError:
            nanoversion = 0
if nanoversion != 0:
    version_nums.append(nanoversion)

if None in version_nums:
    print("Failed to get version number in setup.py.")
    raise

version = ".".join(["%d" % n for n in version_nums[:3]])
if len(version_nums) > 3:
    version += "-%d" % version_nums[3]

extras_require = {
    "doc": [
        "Sphinx==4.5.0",
        "sphinx-autobuild==2021.3.14",
        "sphinxcontrib-bibtex==2.4.2",
        "sphinx-book-theme==0.3.3",
        "myst-parser==0.18.0",
    ],
}

if use_setuptools:
    setup(
        name="spglib",
        version=version,
        cmdclass={"build_ext": CustomBuildExtCommand},
        setup_requires=["numpy", "setuptools>=18.0"],
        license="BSD-3-Clause",
        description="This is the spglib module.",
        long_description=open("README.rst", "rb").read().decode("utf-8"),
        long_description_content_type="text/x-rst",
        author="Atsushi Togo",
        author_email="atz.togo@gmail.com",
        url="http://spglib.github.io/spglib/",
        packages=["spglib"],
        install_requires=[
            "numpy",
        ],
        provides=["spglib"],
        platforms=["all"],
        ext_modules=[extension],
        tests_require=[
            "pyyaml",
        ],
        extras_require=extras_require,
    )
else:
    setup(
        name="spglib",
        version=version,
        license="BSD-3-Clause",
        description="This is the spglib module.",
        long_description=open("README.rst", "rb").read().decode("utf-8"),
        long_description_content_type="text/x-rst",
        author="Atsushi Togo",
        author_email="atz.togo@gmail.com",
        url="http://spglib.github.io/spglib/",
        packages=["spglib"],
        requires=["numpy"],
        provides=["spglib"],
        platforms=["all"],
        ext_modules=[extension],
        tests_require=[
            "pyyaml",
        ],
        extras_require=extras_require,
    )
