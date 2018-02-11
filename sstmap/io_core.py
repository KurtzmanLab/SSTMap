import os
import numpy as np
from string import ascii_uppercase
from io_spatial import do_rotation, rotate_check
from io_helpers import make_grid, are_you_numpy
from collections import OrderedDict
import gzip

### Written by T. Wulsdorf, AG Klebe Marburg University
### 08/2016

class field (object):

    """
    This is a class for operation on generic
    scalar fields. Like GIST objects, they
    must be described in cartesian as well as 
    in fractional space. That means that we need 
    origin, frac2real/real2frac matrix(vector)
    and/or grid spacing vectors vector.

    Tobias Wulsdorf, AG Klebe, 08/2016
    """

    def __init__(self, Bins, Frac2Real=None, Delta=None, Origin=None, Center=None):


        if type(Frac2Real) == type(None) and type(Delta) == type(None):

            raise ValueError("Must provide Frac2Real or Delta.")

        if type(Frac2Real) != type(None) and type(Delta) != type(None):

            raise ValueError("Must provide either Frac2Real or Delta.")

        if type(Frac2Real) == type(None):

            self.delta     = Delta
            self.frac2real = np.eye(3,3) * self.delta

        else:

            self.frac2real = Frac2Real
            self.delta     = np.linalg.norm(self.frac2real, axis=0)

        self.real2frac = np.linalg.inv(self.frac2real)
        self.bins      = Bins

        self.rotation_matrix    = np.eye(3,3)
        self.translation_vector = np.zeros(3)


        if type(Origin) == type(None) and type(Center) == type(None):

            raise ValueError("Must provide origin or Center.")

        if type(Origin) != type(None) and type(Center) != type(None):

            raise ValueError("Must provide either origin or center.")

        if type(Center) == type(None):

            self.origin    = Origin
            self.center    = self.get_real(self.bins/2)

        else:

            self.center    = Center
            #First we need an auxiliary origin at (0,0,0)
            self.origin    = np.zeros(3)
            #Second translate origin according center displacement
            self.origin    = self.center - self.get_real(self.bins/2)

        self.dim = np.array([ np.linalg.norm(self.get_real([self.bins[0], 0., 0.])-self.origin),
                              np.linalg.norm(self.get_real([0., self.bins[1], 0.])-self.origin),
                              np.linalg.norm(self.get_real([0., 0., self.bins[2]])-self.origin)
                              ])


    def translate(self, vector=np.zeros(3)):

        """
        Translatation vector of unit cell origin
        """

        self.translation_vector += vector


    def rotate(self, matrix=np.eye(3,3)):

        """ 
        Rotate the unit cell vectors. 
        """

        rotate_check(matrix)
        self.rotation_matrix = matrix.dot(self.rotation_matrix)


    def translate_global(self, vector=np.zeros(3)):

        """
        Translate global coordinate system
        along vector.
        """

        self.origin += vector


    def rotate_global(self, reference_point=np.zeros(3), matrix=np.eye(3,3)):

        """
        Rotate global coordinate system around
        reference point.
        """

        rotate_check(matrix)
        self.origin = do_rotation(self.origin, reference_point, matrix)

        self.rotate(matrix)

        self.translation_vector = do_rotation(self.translation_vector, np.zeros(3), matrix)


    def get_nice_frac2real(self):

        return self.rotation_matrix.dot(self.frac2real)


    def get_nice_real2frac(self):

        return np.linalg.inv(self.get_nice_frac2real())


    def get_voxel_volume(self):

        """
        Returns the volume per grid voxel.
        """

        return np.absolute(np.cross(self.frac2real[:,0], self.frac2real[:,1]).dot(self.frac2real[:,2]))


    def get_frac(self, real_array):

        #Convert to initial real space by inverse translation and rotation
        initial_reals = do_rotation(real_array, self.origin + self.translation_vector, np.linalg.inv(self.rotation_matrix))

        #Remove origin
        initial_reals -= (self.origin + self.translation_vector)

        #Convert to initial fractional space
        return initial_reals.dot(self.real2frac)


    def get_real(self, frac_array):

        #Convert to real space
        reals = np.array(frac_array).dot(self.frac2real)

        #Perform rotation translation
        return do_rotation(reals, np.zeros(3), self.rotation_matrix) + self.origin + self.translation_vector


    def get_centers(self):

        return self.get_real(make_grid((np.arange(self.bins[0]),\
                                        np.arange(self.bins[1]),\
                                        np.arange(self.bins[2]))))

    def get_centers_real(self):

        return self.get_centers()


    def get_centers_frac(self):

        return make_grid((np.arange(self.bins[0]),\
                          np.arange(self.bins[1]),\
                          np.arange(self.bins[2])))


class gist (field):

  def __init__(self, Bins, Frac2Real=None, Delta=None, Origin=None, Center=None, gist17=False):

    field.__init__(self, Bins, Frac2Real, Delta, Origin, Center)

    if type(gist17) != bool:
      raise IOError("gist17 must be of type bool but is of type %s" %type(gist17))

    ### Compatibility with cpptraj v.17
    self.gist17  = gist17

    ### Number Of nearest neighbours water neighbours
    ### within 3.5 Ang in bulk.
    ### This value is calculted with TIP4PEw
    ### and should be used only in conjunction
    ### with this water model.
    self.bulk_NN = 5.098076
    ### Reference energy interaction energy in kcal/mol.
    ### This value is calculated with TIP4PEW
    self.ref_ene = 11.063656
    ### Reference density in 1/Ang^3
    ### This value is calculated with TIP4PEW as well.
    self.ref_rho = 0.0332

    ### Numpy arrays for storing all quantities
    self.Pop                      = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.gO                       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.gH                       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.dTStrans_dens            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.dTStrans_norm            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.dTSorient_dens           = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.dTSorient_norm           = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.dTSsix_dens              = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.dTSsix_norm              = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Esw_dens                 = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Esw_norm                 = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Eww_dens                 = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Eww_norm_unref           = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Dipole_x_dens            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Dipole_y_dens            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Dipole_z_dens            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Dipole_dens              = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Neighbor_dens            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Neighbor_norm            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self.Order_norm               = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)

    self._tmp_Pop                 = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_gO                  = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_gH                  = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_dTStrans_dens       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_dTStrans_norm       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_dTSorient_dens      = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_dTSorient_norm      = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_dTSsix_dens         = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_dTSsix_norm         = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Esw_norm            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Esw_dens            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Eww_dens            = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Eww_norm_unref      = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Dipole_x_dens       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Dipole_y_dens       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Dipole_z_dens       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Dipole_dens         = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Neighbor_dens       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Neighbor_norm       = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)
    self._tmp_Order_norm          = np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ], dtype=float)

    self._update()


  def _update(self):

    self.Pop            = self._tmp_Pop
    self.gO             = self._tmp_gO
    self.gH             = self._tmp_gH
    self.dTStrans_dens  = self._tmp_dTStrans_dens
    self.dTStrans_norm  = self._tmp_dTStrans_norm
    self.dTSorient_dens = self._tmp_dTSorient_dens
    self.dTSorient_norm = self._tmp_dTSorient_norm
    self.dTSsix_dens    = self._tmp_dTSsix_dens
    self.dTSsix_norm    = self._tmp_dTSsix_norm
    self.Esw_dens       = self._tmp_Esw_dens
    self.Esw_norm       = self._tmp_Esw_norm
    self.Eww_dens       = self._tmp_Eww_dens
    self.Eww_norm_unref = self._tmp_Eww_norm_unref
    self.Dipole_x_dens  = self._tmp_Dipole_x_dens
    self.Dipole_y_dens  = self._tmp_Dipole_y_dens
    self.Dipole_z_dens  = self._tmp_Dipole_z_dens
    self.Dipole_dens    = self._tmp_Dipole_dens
    self.Neighbor_dens  = self._tmp_Neighbor_dens
    self.Neighbor_norm  = self._tmp_Neighbor_norm
    self.Order_norm     = self._tmp_Order_norm


  def cut_round_center(self, bins):

    __doc__="""
    bins is the new self.bins attribute of this class.
    All other attributes (i.e. all scalar fields containing Gist data)
    will be reshaped and transformed such that the position of the
    grid center is preserved. Thereby all edges are cut symmetrically around
    the center.

    Example:
    
    bins=[20,20,20] with self.bins=[50,50,50]
    
    Here, the original self.bins will be reduced to bins. Thereby, the xedge,
    yedge and zedge will be truncated such that [0:15] and [35:49] are cut off.
    """

    if bins.shape[0] != 3:
      print "Target bins array must habe shape (3,)"
    elif self.bins[0] < bins[0] or self.bins[1] < bins[1] or self.bins[2] < bins[2]:
      print "Target bins array must be smaller than original one."
    else:
      cut      = (self.bins - bins)/2
      cut_bins = np.array([[cut[0], cut[0]],
                           [cut[1], cut[1]],
                           [cut[2], cut[2]]])
      self.cut(cut_bins)


  def cut(self, cut_bins):

    __doc__="""
    cut_bins must be an array of shape (3,2), with
    [[lower_cut_x, upper_cut_x],
     [lower_cut_y, upper_cut_y],
     [lower_cut_z, upper_cut_z]]
     marking the number of bins cut out at lower and upper
     boundaries of each dimension.
    """

    self.bins[0] = self.bins[0] - cut_bins[0,0] - cut_bins[0,1]
    self.bins[1] = self.bins[1] - cut_bins[1,0] - cut_bins[1,1]
    self.bins[2] = self.bins[2] - cut_bins[2,0] - cut_bins[2,1]

    self.dim    = self.bins * self.delta
    self.n      = np.copy(self.bins)

    self.origin = np.zeros(3)
    #Second translate origin according to center displacement
    self.origin = self.center - self.get_real(self.bins/2)

    self._tmp_Pop                 = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_gO                  = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_gH                  = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_dTStrans_dens       = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_dTStrans_norm       = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_dTSorient_dens      = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_dTSorient_norm      = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_dTSsix_dens         = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_dTSsix_norm         = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Esw_norm            = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Esw_dens            = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Eww_dens            = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Eww_norm_unref      = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Dipole_x_dens       = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Dipole_y_dens       = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Dipole_z_dens       = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Dipole_dens         = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Neighbor_dens       = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Neighbor_norm       = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)
    self._tmp_Order_norm          = np.array(np.zeros( [ self.bins[0], self.bins[1], self.bins[2] ] ), dtype=float)

    for x_i, x in enumerate(range(cut_bins[0,0], cut_bins[0,0]+self.bins[0])):

      for y_i, y in enumerate(range(cut_bins[1,0], cut_bins[1,0]+self.bins[1])):

        for z_i, z in enumerate(range(cut_bins[2,0], cut_bins[2,0]+self.bins[2])):

          self._tmp_Pop                [x_i][y_i][z_i] = self.Pop            [x,y,z]
          self._tmp_gO                 [x_i][y_i][z_i] = self.gO             [x,y,z]
          self._tmp_gH                 [x_i][y_i][z_i] = self.gH             [x,y,z]
          self._tmp_dTStrans_dens      [x_i][y_i][z_i] = self.dTStrans_dens  [x,y,z]
          self._tmp_dTStrans_norm      [x_i][y_i][z_i] = self.dTStrans_norm  [x,y,z]
          self._tmp_dTSorient_dens     [x_i][y_i][z_i] = self.dTSorient_dens [x,y,z]
          self._tmp_dTSorient_norm     [x_i][y_i][z_i] = self.dTSorient_norm [x,y,z]
          ### In newer vesions of GIST, Six dimensional translational entropy is reported after
          ### dTSorient_norm. For now, we want to skip that.
          if self.gist17:
            self._tmp_dTSsix_dens        [x_i][y_i][z_i] = self.dTSsix_dens    [x,y,z]
            self._tmp_dTSsix_norm        [x_i][y_i][z_i] = self.dTSsix_norm    [x,y,z]
          self._tmp_Esw_dens           [x_i][y_i][z_i] = self.Esw_dens       [x,y,z]
          self._tmp_Esw_norm           [x_i][y_i][z_i] = self.Esw_norm       [x,y,z]
          self._tmp_Eww_dens           [x_i][y_i][z_i] = self.Eww_dens       [x,y,z]
          self._tmp_Eww_norm_unref     [x_i][y_i][z_i] = self.Eww_norm_unref [x,y,z]
          self._tmp_Dipole_x_dens      [x_i][y_i][z_i] = self.Dipole_x_dens  [x,y,z]
          self._tmp_Dipole_y_dens      [x_i][y_i][z_i] = self.Dipole_y_dens  [x,y,z]
          self._tmp_Dipole_z_dens      [x_i][y_i][z_i] = self.Dipole_z_dens  [x,y,z]
          self._tmp_Dipole_dens        [x_i][y_i][z_i] = self.Dipole_dens    [x,y,z]
          self._tmp_Neighbor_dens      [x_i][y_i][z_i] = self.Neighbor_dens  [x,y,z]
          self._tmp_Neighbor_norm      [x_i][y_i][z_i] = self.Neighbor_norm  [x,y,z]
          self._tmp_Order_norm         [x_i][y_i][z_i] = self.Order_norm     [x,y,z]

    self._update()


  def get_nan(self):

    __doc__="""
    Return array that contains True whereever Population array
    is np.nan and False elsewhere.
    """
    tmp = np.zeros(self.bins, dtype=bool)
    tmp[np.isnan(self.Pop)] = True
    
    return tmp


  def get_pop(self):

    __doc__="""
    Return array that containes True whereever Population array
    is greater than zero.
    """

    tmp = np.zeros(self.bins, dtype=bool)
    tmp[np.where(self.Pop > 0)] = True
    
    return tmp


  def write_maps(self, prefix="gist", pymol=True):

    data_dict = OrderedDict()

    data_dict["_Pop.dx"]                = [ self.Pop                                      , 1.0 ]
    data_dict["_gO.dx"]                 = [ self.gO                                       , 4.0 ]
    data_dict["_gH.dx"]                 = [ self.gH                                       , 4.0 ]
    data_dict["_dTStrans_dens.dx"]      = [ self.dTStrans_dens                            , 0.2 ]
    data_dict["_dTStrans_norm.dx"]      = [ self.dTStrans_norm                            , 1.0 ]
    data_dict["_dTSorient_dens.dx"]     = [ self.dTSorient_dens                           , 0.2 ]
    data_dict["_dTSorient_norm.dx"]     = [ self.dTSorient_norm                           , 1.0 ]
    if self.gist17:
      data_dict["_dTSsix_dens.dx"]        = [ self.dTSsix_dens                            , 0.2 ]
      data_dict["_dTSsix_norm.dx"]        = [ self.dTSsix_norm                            , 1.0 ]
    data_dict["_Esw_dens.dx"]           = [ self.Esw_dens                                 , 0.2 ]
    data_dict["_Esw_norm.dx"]           = [ self.Esw_norm                                 , 1.0 ]
    data_dict["_Eww_dens.dx"]           = [ self.Eww_dens                                 , 0.2 ]
    data_dict["_Eww_norm_unref.dx"]     = [ self.Eww_norm_unref                           , 1.0 ]
    data_dict["_Eww_norm_ref.dx"]       = [ self.ref_ene-self.Eww_norm_unref              , 1.0 ]
    data_dict["_Eww_norm_ref_dens.dx"]  = [ (self.ref_ene-self.Eww_norm_unref)*\
                                             self.gO*self.ref_rho                         , 1.0 ]
    data_dict["_Dipole_x_dens.dx"]      = [ self.Dipole_x_dens                            , 1.0 ]
    data_dict["_Dipole_y_dens.dx"]      = [ self.Dipole_y_dens                            , 1.0 ]
    data_dict["_Dipole_z_dens.dx"]      = [ self.Dipole_z_dens                            , 1.0 ]
    data_dict["_Dipole_dens.dx"]        = [ self.Dipole_dens                              , 1.0 ]
    data_dict["_Neighbor_dens.dx"]      = [ self.Neighbor_dens                            , 0.5 ]
    data_dict["_Neighbor_norm.dx"]      = [ self.Neighbor_norm                            , 1.0 ]
    data_dict["_Order_norm.dx"]         = [ self.Order_norm                               , 5.5 ]
    data_dict["_dTS_dens.dx"]           = [ self.dTStrans_dens+self.dTSorient_dens        , 0.2 ]
    data_dict["_dTS_norm.dx"]           = [ self.dTStrans_norm+self.dTSorient_norm        , 1.0 ]
    data_dict["_E_dens.dx"]             = [ self.Esw_dens+self.Eww_dens                   , 0.2 ]
    data_dict["_E_norm.dx"]             = [ self.Esw_norm+self.Eww_norm_unref             , 1.0 ]
    data_dict["_Neighbor_loss_norm.dx"] = [ self.Neighbor_norm - self.bulk_NN             , 0.5 ]

    if pymol:

      pymol_string = ""
      pymol_string += "from pymol import cmd\n"
      pymol_string += "from collections import OrderedDict\n"
      pymol_string += "\n"

    for name, data in data_dict.items():

      write_files(Frac2Real=self.get_nice_frac2real(), Bins=self.bins, Origin=self.origin, Value=data[0], Format="DX", Filename=prefix+name, Nan_fill=-999)

      if pymol:

        new_name      = str(prefix+name).replace(".dx", "")

        pymol_string += "### %s ###\n" %new_name
        pymol_string += "cmd.load(\"./%s\")\n" %(prefix+name)
        pymol_string += "cmd.isomesh(\"%s\", \"%s\", level=%s)\n" %(new_name+"_map", new_name, data[1])
        pymol_string += "cmd.map_double(\"%s\")\n" %(new_name)
        pymol_string += "\n"

    if pymol:

      pymol_string += "cmd.disable(\"*_map\")\n"
      pymol_string += "cmd.do(\"color blue, *_map\")\n"
      pymol_string += "cmd.do(\"set mesh_negative_color, red\")\n"
      pymol_string += "cmd.do(\"set mesh_negative_visible\")\n"

      with open(prefix+"_pymol.py", "w") as f:

        f.write(pymol_string)


class loadgist(gist):

  def __init__(self, Path, gist17=False):

    gist.__init__(self, Bins=np.array([50,50,50]), Origin=np.array([0,0,0]), Delta=np.array([0.5,0.5,0.5]), gist17=gist17)

    if not os.path.exists(Path):
      raise IOError("File %s not found." %Path)

    self.path    = Path

    if Path.endswith(".gz"):
      map_file_ref = gzip.open(Path,"r")
    else:
      map_file_ref = open(Path,"r")
    map_file     = map_file_ref.readlines()
    map_file_ref.close()

    ##### Read in the gist-out file produced by the
    ##### the gist functionality of cpptraj
 
    start_row = -1

    for i, item in enumerate(map_file):

      if len(item.rstrip().split()) > 1 and item.rstrip().split()[0] == 'voxel' \
                                        and item.rstrip().split()[1] == 'xcoord':

        start_row = i + 1
        break

    z_start = map_file[start_row].rstrip().split()[3]
    y_start = map_file[start_row].rstrip().split()[2]
    x_start = map_file[start_row].rstrip().split()[1]

    found_bins_z = False
    found_bins_y = False
    found_bins_x = False


    # Find out grid dimensions and bin spacing
    # We assume that the cell is completely rectangular
    #
    # The coordinate data is ordered in opendx like fashion:
    # The z coordinate is running fastest, then y coordinate,
    # then x coordinate. E.g.:
    # (x_0, y_0, z_0), (x_0, y_0, z_1), (x_0, y_0, z_2), ...

    for i, line in enumerate(map_file[start_row+1:]):

      if found_bins_z and not found_bins_y and line.rstrip().split()[2] == y_start:

        self.bins[1] = (i+1) / self.bins[2]

        break

      if not found_bins_z and line.rstrip().split()[3] == z_start:

        self.bins[2] = i + 1

        found_bins_z = True

    ### Now we fill the grids...

    self.bins[0] = (len(map_file) - start_row) / (self.bins[2] * self.bins[1])

    # self.n is only here for historical reasons.
    # ... and for compatibility with xplor maps!
    self.n = np.copy(self.bins)

    self.delta  = np.array( [ float(map_file[ start_row + 1 + self.bins[0] * self.bins[1]].rstrip().split()[1] )\
                                                            - float(map_file[ start_row ].rstrip().split()[1] ),
                              float(map_file[ start_row + 1 + self.bins[0] ].rstrip().split()[2] )\
                                                            - float(map_file[ start_row ].rstrip().split()[2] ),
                              float(map_file[ start_row + 1 ].rstrip().split()[3] )\
                                                            - float(map_file[ start_row ].rstrip().split()[3] ) ] )

    self.dim    = self.bins * self.delta

    self.origin = -self.delta * 0.5 + np.array( [ float(map_file[ start_row ].rstrip().split()[1] ),
                                                  float(map_file[ start_row ].rstrip().split()[2] ),
                                                  float(map_file[ start_row ].rstrip().split()[3] ) ] )

    ### We assume that all cell angles are 90 deg.
    ### Cell edges can be of different length
    self.frac2real = np.array( [ [ self.dim[0]  / self.n[0], 0.0,                      0.0                       ],
                                 [ 0.0,                      self.dim[1]  / self.n[1], 0.0                       ],
                                 [ 0.0,                      0.0,                      self.dim[2]  / self.n[2]  ] 
                                ] 
                            )

    self.real2frac          = np.linalg.inv(self.frac2real)
    self.rotation_matrix    = np.eye(3,3)
    self.translation_vector = np.zeros(3)
    self.center             = self.get_real(self.bins/2)

    ###### now we read all the gist-specific data...

    i = 0

    if self.gist17:
      skip_col=2
    else:
      skip_col=0

    for x in range(0,self.bins[0]):

      for y in range(0,self.bins[1]):

        for z in range(0,self.bins[2]):

          self._tmp_Pop                [x][y][z] = float(map_file[start_row + i].rstrip().split()[4])
          self._tmp_gO                 [x][y][z] = float(map_file[start_row + i].rstrip().split()[5])
          self._tmp_gH                 [x][y][z] = float(map_file[start_row + i].rstrip().split()[6])
          self._tmp_dTStrans_dens      [x][y][z] = float(map_file[start_row + i].rstrip().split()[7])
          self._tmp_dTStrans_norm      [x][y][z] = float(map_file[start_row + i].rstrip().split()[8])
          self._tmp_dTSorient_dens     [x][y][z] = float(map_file[start_row + i].rstrip().split()[9])
          self._tmp_dTSorient_norm     [x][y][z] = float(map_file[start_row + i].rstrip().split()[10])
          ### In newer vesions of GIST, Six dimensional translational entropy is reported after
          ### dTSorient_norm. For now, we want to skip that.
          if self.gist17:
            self._tmp_dTSsix_dens        [x][y][z] = float(map_file[start_row + i].rstrip().split()[11])
            self._tmp_dTSsix_norm        [x][y][z] = float(map_file[start_row + i].rstrip().split()[12])
          self._tmp_Esw_dens           [x][y][z] = float(map_file[start_row + i].rstrip().split()[11+skip_col])
          self._tmp_Esw_norm           [x][y][z] = float(map_file[start_row + i].rstrip().split()[12+skip_col])
          self._tmp_Eww_dens           [x][y][z] = float(map_file[start_row + i].rstrip().split()[13+skip_col])
          self._tmp_Eww_norm_unref     [x][y][z] = float(map_file[start_row + i].rstrip().split()[14+skip_col])
          self._tmp_Dipole_x_dens      [x][y][z] = float(map_file[start_row + i].rstrip().split()[15+skip_col])
          self._tmp_Dipole_y_dens      [x][y][z] = float(map_file[start_row + i].rstrip().split()[16+skip_col])
          self._tmp_Dipole_z_dens      [x][y][z] = float(map_file[start_row + i].rstrip().split()[17+skip_col])
          self._tmp_Dipole_dens        [x][y][z] = float(map_file[start_row + i].rstrip().split()[18+skip_col])
          self._tmp_Neighbor_dens      [x][y][z] = float(map_file[start_row + i].rstrip().split()[19+skip_col])
          self._tmp_Neighbor_norm      [x][y][z] = float(map_file[start_row + i].rstrip().split()[20+skip_col])
          self._tmp_Order_norm         [x][y][z] = float(map_file[start_row + i].rstrip().split()[21+skip_col])

          i += 1

    self._update()

class write_files (object):

    def __init__(self, Delta=None, Frac2Real=None, Bins=None, Origin=None, \
                 Value=None, XYZ=None, X=None, Y=None, Z=None, Format='PDB', \
                 Filename=None, Nan_fill=-1.0):

        """
        This class can write different file types.
        currently only dx and pdb are supported.
        """

        self._delta     = Delta
        self._frac2real = Frac2Real
        self._bins      = Bins
        self._origin    = Origin
        self._value     = Value
        self._x         = X
        self._y         = Y
        self._z         = Z
        self._format    = Format
        self._filename  = Filename
        self._xyz       = XYZ
        self._nan_fill  = Nan_fill

        if type(self._filename) != str:

            self._filename  = 'output.'
            self._filename  += self._format

        self._writers = {
                        'PDB'  : self._write_PDB,
                        'DX'   : self._write_DX,
                        'GIST' : self._write_GIST
                        }

        data = self._writers[self._format]()

        o = open(self._filename, "w")
        o.write(data)
        o.close()

    def _merge_x_y_z(self):

        return np.stack( ( self._x, self._y, self._z ), axis=1 )


    def _write_PDB(self):

        """
        Write a PDB file.
        This is intended for debugging. It writes all atoms
        as HETATM of element X with resname MAP.
        """

        if are_you_numpy(self._xyz):

            if self._xyz.shape[-1] != 3:

                raise TypeError(
                    "XYZ array has wrong shape.")

        else:

            if not ( are_you_numpy(self._x) or are_you_numpy(self._y) or are_you_numpy(self._z) ):

                raise TypeError(
                    "If XYZ is not given, x,y and z coordinates\
                     must be given in separate arrays.")

            else:

                self._xyz = self._merge_x_y_z()

        if self._value == None:

            self._value = np.zeros( len(self._xyz), dtype=float )

        data = 'REMARK File written by write_files.py\n'

        for xyz_i, xyz in enumerate(self._xyz):

            #iterate over uppercase letters
            chain_id    = ascii_uppercase[( len(str(xyz_i+1)) / 5 )]

            atom_counts = xyz_i - ( len(str(xyz_i+1)) / 6 ) * 100000
            resi_counts = xyz_i - ( len(str(xyz_i+1)) / 5 ) * 10000
            data += \
            '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          \n' \
            %('HETATM',atom_counts+1,'X','', 'MAP', chain_id, resi_counts+1, '', xyz[0], xyz[1], xyz[2], 0.00, float( self._value[xyz_i] ) )

        data += 'END\n'

        return data


    def _write_DX(self):

        """
        Writes DX files according to openDX standard.
        """

        if not ( are_you_numpy(self._origin) or are_you_numpy(self._bins) ):

            raise TypeError(
            "Origin and bins must be given.")

        #This means not (a XOR b) or not (a or b)
        if are_you_numpy(self._delta) == are_you_numpy(self._frac2real) :

            raise TypeError(
            "Either delta or frac2real must be given.")

        if are_you_numpy(self._delta):

            self._frac2real = np.zeros((3,3), dtype=float)

            np.fill_diagonal(self._frac2real, self._delta)

        data = '''object 1 class gridpositions counts %d %d %d
origin %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
object 2 class gridconnections counts %d %d %d
object 3 class array type float rank 0 items %d data follows
''' %(self._bins[0], self._bins[1], self._bins[2],\
      self._origin[0], self._origin[1], self._origin[2],\
      self._frac2real[0][0], self._frac2real[0][1], self._frac2real[0][2],\
      self._frac2real[1][0], self._frac2real[1][1], self._frac2real[1][2],\
      self._frac2real[2][0], self._frac2real[2][1], self._frac2real[2][2],\
      self._bins[0],   self._bins[1],   self._bins[2],\
      self._bins[2] * self._bins[1] * self._bins[0])

        i = 0
        for x_i in range(0, self._bins[0]):

            for y_i in range(0, self._bins[1]):

                for z_i in range(0, self._bins[2]):

                    ### writing an integer instead of float
                    ### saves us some disk space
                    if np.isnan(self._value[x_i][y_i][z_i]):

                        data += str(self._nan_fill) + " "

                    else:

                        if self._value[x_i][y_i][z_i] == 0.0:

                            data += "0 " 

                        else:

                            data += str(self._value[x_i][y_i][z_i]) + ' '
                
                    i += 1

                    if i == 3:

                        data += '\n'
                        i = 0
        return data

    def _write_GIST(self):

        """
        To be implemented...
        """

        pass

class PDB(object):

  """
  Class that reads a pdb file and provides pdb type data structure.
  """

  def __init__(self, Path):

    self.path = Path

    self.crd  = list()
    self.B    = list()

    with open(self.path, "r") as PDB_file:

      for i, line in enumerate(PDB_file):

        if not (line[0:6].rstrip() == 'ATOM' or line[0:6].rstrip() == 'HETATM'):

          continue

        if i <= 9999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[30:38]))
          self.crd[-1].append(float(line.rstrip()[38:46]))
          self.crd[-1].append(float(line.rstrip()[46:54]))

          #B-Factors
          self.B.append(line.rstrip()[54:59])

        if 9999 < i <= 99999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[31:39]))
          self.crd[-1].append(float(line.rstrip()[39:47]))
          self.crd[-1].append(float(line.rstrip()[47:55]))

          #B-Factors
          self.B.append(line.rstrip()[55:60])

        if i > 99999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[33:41]))
          self.crd[-1].append(float(line.rstrip()[41:49]))
          self.crd[-1].append(float(line.rstrip()[49:57]))

          #B-Factors
          self.B.append(line.rstrip()[57:62])

    self.crd  = np.array(self.crd)
    self.B    = np.array(self.B)

def guess_field(crds, delta=np.array([0.5,0.5,0.5])):

    _c = np.mean(crds, axis=0)

    _min = np.min((_c - crds), axis=0)
    _max = np.max((_c - crds), axis=0)

    _b = np.rint(np.abs(_max - _min)/delta + (5.0 / delta) )

    del _min, _max

    return field(Bins=_b, Delta=delta, Center=_c)