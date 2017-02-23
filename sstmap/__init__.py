#############################################################################
# SSTMap: A Python library for the calculation of water structure and 
#         thermodynamics on solute surfaces from molecular dynamics 
#         trajectories.
# Copyright 2016-2017 Lehman College City University of New York and the Authors
#
# Authors: Kamran Haider
# Contributors: Steven Ramsay, Anthony Cruz Balberdy
#
# SSTMap is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################
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

__author__ = "Kamran Haider"
__license__ = "LGPL 2.1"
__maintainer__ = "Kamran Haider"
__email__ = "kamranhaider.mb@gmail.com"


from sstmap import site_water_analysis, grid_water_analysis, utils
from sstmap.site_water_analysis import SiteWaterAnalysis
from sstmap.grid_water_analysis import GridWaterAnalysis
