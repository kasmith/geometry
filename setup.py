#!/usr/bin/env python
import distutils
from distutils.core import setup

setup(name = "geometry",
      version = '0.5',
      description = "Tools for complex geometric parsing of scenes and shapes",
      author = "Kevin A Smith",
      author_email= "k2smith@mit.edu",
      url = "Https://github.com/kasmith/geometry",
      packages = ["geometry", "geometry.helpers"],
      requires = ["numpy"])
