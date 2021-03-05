import sys, os, string, numba
import numpy as np
import pygal
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, thermal, border, magnitude_1
from plaq_v_cfg import Nstart, Nend
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import gauge_latticeqcd
import matplotlib as cm
import tools_v1
import Relative_tools