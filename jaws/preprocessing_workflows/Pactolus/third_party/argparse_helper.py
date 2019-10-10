"""
COPYRIGHT: This module has been extracted from BASTet (Berkeley Analysis and Storage Toolkit).
           In BASTet the sources are in omsi/workflow/common.py and in part of omsi/datastructures/analysis_data.
           The sources have been modified to remove dependencies on other parts of BASTet

*** Copyright Notice ***

BASTet  Copyright (c) 2015, The Regents of the University of California, through Lawrence Berkeley National Laboratory
(subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's
Innovation & Partnerships Office at IPO@lbl.gov.

NOTICE.  This software was developed under funding from the U.S. Department of Energy.  As such, the U.S. Government
has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license
in the Software to reproduce, prepare derivative works, and perform publicly and display publicly.  Beginning five (5)
years after the date permission to assert copyright is obtained from the U.S. Department of Energy, and subject to
any subsequent five (5) year renewals, the U.S. Government is granted for itself and others acting on its behalf a
paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works,
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

*** License***

Lawrence Berkeley National Laboratory

ACADEMIC, NON-PROFIT, NON-COMMERCIAL USE ONLY, LICENSE

BASTet Copyright (c) 2015, The Regents of the University of California, through Lawrence Berkeley National Laboratory
(subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this list of conditions and the
    following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
    following disclaimer in the documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy
    nor the names of its contributors may be used to endorse or promote products derived from this software without
    specific prior written permission.

(4) Use of the software, in source or binary form is for academic, non-profit, NON-COMMERCIAL USE, purposes ONLY.
    All commercial use rights for the software are hereby reserved. A separate commercial use license is available
    from Lawrence Berkeley National Laboratory.

(5) In the event you create any bug fixes, patches, upgrades, updates, modifications, derivative works or enhancements
    to the source code or binary code of the software ("Enhancements") you hereby grant The Regents of the University
    of California and the U.S. Government a paid-up, non-exclusive, irrevocable, worldwide license in the Enhancements
    to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly,
    and to permit others to do so.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

#########################################################
#      omsi/workflow/common.py                          #
#########################################################

import argparse
import numpy as np
import ast

class RawDescriptionDefaultHelpArgParseFormatter(argparse.ArgumentDefaultsHelpFormatter,
                                                 argparse.RawDescriptionHelpFormatter):
    """
    Simple derived formatter class for use with argparse
    """
    pass


#########################################################
#      Part of omsi/analysis/analysis_data.py           #
#########################################################
class data_dtypes(dict):
    """
    Class specifying basic function for specifying common
    data types used as part of an analysis.
    """
    @staticmethod
    def get_dtypes():
        """
        Get a list of available data type specifications
        """
        dtypes = {'int': int,
                  'float': float,
                  'long': long,
                  'complex': complex,
                  'bool': data_dtypes.bool_type,
                  'str': str,
                  'unicode': unicode,
                  'ndarray': data_dtypes.ndarray}
        return dtypes

    @staticmethod
    def bool_type(argument):
        """
        Implement conversion of boolean input parameters since
        arparse (or bool, depending on the point of view), do not
        handle bool as a type in an intuitive fashion.

        :param argument: The argument to be parsed to a boolean
        :return: The converted value
        """
        try:
            bool(int(argument))
        except ValueError:
            if argument in ('TRUE', 'true', 'True', 't', 'T'):
                return True
            elif argument in ('FALSE', 'false', 'False', 'f', 'F'):
                return False
            else:
                raise ValueError('Parameter could not be converted to type bool')

    @staticmethod
    def ndarray(argument):
        """
        This dtype may be used to indicate numpy ndarrays as
        well as h5py arrays or omsi_dependencies

        :param argument: The argument to be parsed to ndarray

        :return: The converted ndarray
        """
        # from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        # from omsi.dataformat.omsi_file.msidata import omsi_file_msidata
        # from omsi.dataformat.omsi_file.common import omsi_file_common
        if isinstance(argument, basestring):
            try:
                return np.asarray(ast.literal_eval(argument))
            except (ValueError, SyntaxError):
                # omsi_out_object = omsi_file_common.get_omsi_object(h5py_object=argument)
                # if omsi_out_object is not None:
                #     return omsi_out_object
                # else:
                #
                raise ValueError('String could not be converted to valid ndarray. This may be ' +
                                 'due to, e.g., a syntax error or the file may not exists')
        # elif isinstance(argument, dependency_dict) or \
        #        isinstance(argument, h5py.Dataset) or isinstance(argument, h5py.Group) or \
        #        isinstance(argument, omsi_file_analysis) or \
        #        isinstance(argument, omsi_file_msidata):
        #    return argument
        elif argument is None:
            return None
        return np.asarray(argument)
