##############################################################################
#  SSTMap: A Python library for the calculation of water structure and
#         thermodynamics on solute surfaces from molecular dynamics
#         trajectories.
# MIT License
# Copyright 2016-2017 Lehman College City University of New York and the Authors
#
# Authors: Kamran Haider, Steven Ramsey, Anthony Cruz Balberdy, Tom Kurtzman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################
"""
This module contains implementation of a parent class for water analysis in 
molecular dynamics simulation trajectories. This class provides methods for 
index all atoms in the simulation, calculations of the energy and hydrogen 
bonding of water molecules with other atoms in the system.

Please reference the following if you use this code in your research:
[1] Haider K, Wickstrom L, Ramsey S, Gilson MK and Kurtzman T. Enthalpic Breakdown 
of Water Structure on Protein Active Site Surfaces. J Phys Chem B. 120:8743-8756, 
2016. http://dx.doi.org/10.1021/acs.jpcb.6b01094.
"""

__author__ = "Kamran Haider, Anthony Cruz, Steven Ramsey, Tobias Wulsdorf, Tom Kurtzman"
__license__ = "MIT"
__maintainer__ = "Kamran Haider"
__email__ = "kamranhaider.mb@gmail.com"


from sstmap import site_water_analysis, grid_water_analysis, utils
from sstmap.site_water_analysis import SiteWaterAnalysis
from sstmap.grid_water_analysis import GridWaterAnalysis
