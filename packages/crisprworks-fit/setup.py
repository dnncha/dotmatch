from setuptools import Extension, setup

setup(ext_modules=[Extension(
    "crisprworks_fit._native", ["src/crisprworks_fit/_native.c"],
    extra_compile_args=["-O3", "-fno-fast-math", "-ffp-contract=off"], optional=True,
)])
