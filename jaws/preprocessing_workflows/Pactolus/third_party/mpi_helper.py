"""
Module used to ease the use of MPI and distributed parallel implementations using MPI
"""

"""
COPYRIGHT: This module has been extracted from BASTet (Berkeley Analysis and Storage Toolkit). (omsi/shared/mpi_helper.py)

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

try:
    def test_mpi_available():
        """
        This function import MPI in a seperate process to safely check if
        MPI is available. This precaution is necessary as on Cray systems
        importing MPI can lead to a crash on, e.g., login nodes where the
        use of MPI is not permitted. By executing the import in a separate
        process we avoid crashing the main process and we can safely check
        whether the process aborted or not.

        :return: False if the import failed, otherwise return True
        """
        import sys
        from subprocess import Popen, PIPE
        process = Popen('%s -c "from mpi4py import MPI as mpi"'%( \
                         sys.executable), shell=True, stderr=PIPE, stdout=PIPE)
        import_failed = process.wait()
        return not import_failed

    MPI_AVAILABLE = test_mpi_available()

    if MPI_AVAILABLE:
        from mpi4py import MPI
except ImportError:
    MPI_AVAILABLE = False

if not MPI_AVAILABLE:
    try:
        from omsi.shared.log import log_helper
        log_helper.warning(__name__, "MPI not available. Running in serial.")
    except:
        print "MPI not available. Running in serial."

import numpy as np
import itertools
import warnings
import time


class parallel_over_axes(object):
    """
    Helper class used to parallelize the execution of a function using MPI by splitting the
    input data into sub-blocks along a given set of axes.

    :ivar task_function: The function we should run.
    :ivar task_function_params: Dict with the input parameters for the function.
        may be None or {} if no parameters are needed.
    :ivar main_data: Dataset over which we should parallelize
    :ivar split_axes:  List of integer axis indicies over which we should parallelize
    :ivar main_data_param_name: The name of data input parameter of the task function
    :ivar root: The master MPI rank (Default=0)
    :ivar schedule: The task scheduling schema to be used (see parallel_over_axes.SCHEDULES
    :ivar collect_output: Should we collect all the output from the ranks on the master rank?
    :ivar schedule: The parallelization schedule to be used. See also parallel_over_axes.schedule
    :ivar result: The result form the task_function. If self.__data_collected  is set and we are the root
        then this will a list with the the output of all tasks
    :ivar blocks: List with tuples describing the selected subset of data processed by the given block task.
        If self.__data_collected is set and we are the root rank then this is a list of all the blocks
        processed by each rank.
    :ivar block_times: List of times in seconds used to process the data block with the given index.
        NOTE: The block times include also any required communications and other operations to initialize
        and complete the task, and not just the execution of the task function itself.
    :ivar run_time: Float time in seconds for executing the run function.
    :ivar comm: The MPI communicator used for the parallelization. Default value is MPI.COMM_WORLD

    """
    SCHEDULES = {'STATIC_1D': 'STATIC_1D',
                 'STATIC': 'STATIC',
                 'DYNAMIC': 'DYNAMIC'}

    MPI_MESSAGE_TAGS = {'RANK_MSG': 11,
                        'BLOCK_MSG': 12,
                        'COLLECT_MSG': 13}

    def __init__(self,
                 task_function,
                 task_function_params,
                 main_data,
                 split_axes,
                 main_data_param_name,
                 schedule=SCHEDULES['STATIC_1D'],
                 root=0,
                 comm=None):
        """

        :param task_function: The function we should run.
        :param task_function_params: Dict with the input parameters for the function.
            May be None or {} if no parameters are needed.
        :param main_data: Dataset over which we should parallelize
        :param split_axes:  List of integer axis indicies over which we should parallelize
        :param main_data_param_name: The name of data input parameter of the task function
        :param root: The master MPI rank (Default=0)
        :param schedule: The task scheduling schema to be used (see parallel_over_axes.SCHEDULES
        :param comm: The MPI communicator used for the parallelization. Default value is None, in which case
            MPI.COMM_WORLD is used

        """
        if not is_mpi_available():
            raise ValueError("MPI is not available. MPI is required for parallel execution.")
        self.task_function = task_function
        self.schedule = schedule
        self.split_axes = split_axes
        if isinstance(self.split_axes, int):  # Make sure that split-axis is a list not just a single index
            self.split_axes = [self.split_axes, ]
        self.main_data = main_data
        self.main_data_param_name = main_data_param_name
        self.task_function_params = task_function_params
        if self.task_function_params is None:
            self.task_function_params = {}
        self.root = root
        self.result = None
        self.blocks = None
        self.block_times = None
        self.run_time = None
        self.__data_collected = False
        self.comm = get_comm_world() if comm is None else comm

    def run(self):
        """
        Call this function to run the function in parallel.

        :return: Tuple with the following elements:

            1) List with the results from the local execution of the task_function. Each
               entry is the result from one return of the task_function. In the case of static
               execution, this is always a list of length 1.
            2) List of block_indexes. Each block_index is a tuple with the selection used to
               divide the data into sub-blocks. In the case of static decomposition we have
               a range slice object along the axes used for decomposition whereas in the
               case of dynamic scheduling we usually have single integer point selections
               for each task.

        """
        try:
            from omsi.shared.log import log_helper
        except ImportError:
            from pactolus.third_party.log import log_helper
        start_time = time.time()
        self.__data_collected = False
        if self.schedule == self.SCHEDULES['DYNAMIC']:
            result = self.__run_dynamic()
        elif self.schedule == self.SCHEDULES['STATIC_1D'] or self.schedule == self.SCHEDULES['STATIC']:
            result = self.__run_static_1D()
        else:
            log_helper.error(__name__, "Invalid scheduling scheme given: " + str(self.schedule))
            raise ValueError("Invalid scheduling scheme given: " + str(self.schedule))
        end_time = time.time()
        self.run_time = end_time - start_time
        return result

    def collect_data(self, force_collect=False):
        """
        Collect the results from the parallel execution to the self.root rank.

        NOTE: On the root the self.result, self.blocks, and self.block_times variables are
              updated with the collected data as well and self.__data_collected will be set

        NOTE: If the data has already been collected previously (ie., collect_data has been called
            before), then the collection will not be performed again, unless force_collect is set.

        :param force_collect: Set this parameter to force that data collection is performed again.
            By default the collect_data is performed only once for each time the run(..) function
            is called and the results are reused to ensure consistent data structures. We can
            force that collect will be reexecuted anyways by setting force_collect.

        :return: On worker ranks (i.e., MPI_RANK!=self.root) this is simply the
            self.result and self.blocks containing the result created by run function.
            On the root rank (i.e., MPI_RANK!=self.root) this is a tuple of two lists
            containing the combined data of all  self.result and self.blocks from all ranks respectively.

        """
        try:
            from omsi.shared.log import log_helper
        except ImportError:
            from pactolus.third_party.log import log_helper
        # If we have collected the data already then we don't need to do it again
        if self.__data_collected and not force_collect:
            return self.result, self.blocks

        # Collect the output
        rank = get_rank(comm=self.comm)
        start_time = time.time()
        if rank == self.root:
            log_helper.info(__name__, "COLLECTING RESULTS")
        # Collect the data, blocks, and block_times from all ranks
        collected_data = self.comm.gather(self.result, root=self.root)
        collected_blocks = self.comm.gather(self.blocks, root=self.root)
        # Save the data to self.result, self.block, self.block_times if we are the root
        if rank == self.root:
            # Merge the results from all the processes into a single result and blocks list
            # rather than having a list of lists of results
            self.result = list(itertools.chain.from_iterable(collected_data))
            self.blocks = list(itertools.chain.from_iterable(collected_blocks))

        # Record the time we used to collect the data
        end_time = time.time()
        run_time = end_time - start_time
        if rank == self.root:
            log_helper.info(__name__, "TIME FOR COLLECTING DATA FROM ALL TASKS: " + str(run_time))
        # Return the result
        self.__data_collected = True
        return self.result, self.blocks

    def __run_static_1D(self):
        """
        Run the task function using a static task decomposition schema.

        The data is divided into sub-blocks along the largest split_axis

        :return: Tuple with the following elements:

            1) List with the results from the local execution of the task_function. Each
               entry is the result from one return of the task_function. In the case of static
               execution, this is always a list of length 1.
            2) List of block_indexes. Each block_index is a tuple with the selection used to
               divide the data into sub-blocks. In the case of static decomposition we have
               a range slice object along the axes used for decomposition.

        """
        try:
            from omsi.shared.log import log_helper
        except ImportError:
            from pactolus.third_party.log import log_helper
        start_time = time.time()
        # Get MPI parameters
        rank = get_rank(comm=self.comm)
        size = get_size(comm=self.comm)

        # Get data shape parameters and compute the data blocks
        # Determine the longest axis along which we can split the data
        axes_shapes = np.asarray(self.main_data.shape)[self.split_axes]
        total_num_subblocks = np.prod(axes_shapes)
        if total_num_subblocks < size:
            size = total_num_subblocks
            if rank == self.root:
                log_helper.info(__name__,
                                "Insufficient number of blocks for number of MPI ranks. Some ranks will remain idle")
        axes_sort_index = np.argsort(axes_shapes)[::-1]
        split_axis = self.split_axes[axes_sort_index[0]]
        split_axis_size = axes_shapes[split_axis]
        if split_axis_size < size:
            raise NotImplementedError("STATIC scheduling currently parallelizes only over one axis, " +
                                      "and the largest axis is too small to fill all MPI tasks")
        # Determine the size of 1D block
        block_size = int(split_axis_size / float(size) + 0.5)
        if block_size * size > split_axis_size and block_size > 1:
            block_size -= 1

        # Compute a block for every rank
        self.blocks = [slice(None)] * len(self.main_data.shape)
        start_index = rank * block_size
        stop_index = start_index + block_size
        if rank == (size-1):
            if stop_index != split_axis_size:
                stop_index = split_axis_size
        self.blocks[axes_sort_index[0]] = slice(start_index, stop_index)
        self.blocks = tuple(self.blocks)
        log_helper.info(__name__, "Rank: " + str(rank) + " Block: " + str(self.blocks))

        # Execute the task_function on the given data block
        task_params = self.task_function_params
        task_params[self.main_data_param_name] = self.main_data[self.blocks]
        self.result = self.task_function(**task_params)

        end_time = time.time()
        run_time = end_time - start_time
        self.block_times = [run_time, ]
        log_helper.info(__name__, "TIME FOR PROCESSING THE DATA BLOCK: " + str(run_time))

        # Return the output
        self.result = [self.result, ]
        self.blocks = [self.blocks, ]
        return self.result, self.blocks

    def __run_dynamic(self):
        """
        Run the task function using dynamic task scheduling.

        The root rank divides the data into sub-tasks and sends the tasks to available MPI
        processes on request.

        :return: Tuple with the following elements:

            1) List with the results from the local execution of the task_function. Each
               entry is the result from one return of the task_function.
            2) List of block_indexes. Each block_index is a tuple with the selection used to
               divide the data into sub-blocks. In the case of static decomposition we have
               a range slice object along the axes used for decomposition.

        """
        try:
            from omsi.shared.log import log_helper
        except ImportError:
            from pactolus.third_party.log import log_helper
        import time
        rank = get_rank(comm=self.comm)
        size = get_size(comm=self.comm)

        if size < 2:
            warnings.warn('DYNAMIC task scheduling requires at least 2 MPI ranks. Using STATIC scheduling instead.')
            return self.__run_static_1D()

        # We are the controlling rank
        if rank == self.root:
            self.result = []
            self.blocks = []
            self.block_times = []
            # Get data shape parameters and compute the data blocks
            axes_shapes = np.asarray(self.main_data.shape)[self.split_axes]
            total_num_subblocks = np.prod(axes_shapes)
            if total_num_subblocks < size:
                if rank == self.root:
                    warnings.warn("Insufficient number of blocks for number of MPI ranks. Some ranks will remain idle")

            # Compute the list of all possible blocks
            base_blocks = [[slice(None)]] * len(self.main_data.shape)
            for axis_index in self.split_axes:
                base_blocks[axis_index] = range(self.main_data.shape[axis_index])
            block_tuples = itertools.product(*base_blocks)

            # Communicate blocks with task ranks
            log_helper.info(__name__, "PROCESSING DATA BLOCKS")
            start_time = time.time()
            block_index = 0
            for block_selection in block_tuples:
                request_rank = self.comm.recv(source=MPI.ANY_SOURCE, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                self.comm.send((block_index, block_selection),
                               dest=request_rank,
                               tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                block_index += 1
                if (block_index % 100) == 0:
                    log_helper.debug(__name__, str((block_index, total_num_subblocks, request_rank)))
            end_time = time.time()
            run_time = end_time - start_time
            log_helper.info(__name__, "TIME FOR SCHEDULING ALL TASKS: " + str(run_time))
            start_time = time.time()
            log_helper.info(__name__, "FINALIZING")
            # Terminate all ranks and receive all data from the different ranks if requested
            all_ranks_status = np.zeros(size, 'bool')
            all_ranks_status[self.root] = True
            while not np.all(all_ranks_status):

                request_rank = self.comm.recv(source=MPI.ANY_SOURCE, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                self.comm.send((None, None), dest=request_rank, tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                all_ranks_status[request_rank] = True

            end_time = time.time()
            run_time = end_time - start_time
            log_helper.info(__name__, "TIME FOR FINALIZING TASKS: " + str(run_time))

        # We are a rank that has to run tasks
        else:
            # Request a new data block
            self.result = []
            self.blocks = []
            self.block_times = []
            while True:
                start_time = time.time()
                self.comm.send(rank, dest=self.root, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                block_index, block_selection = self.comm.recv(source=self.root, tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                if block_index is None:
                    break
                # Execute the task_function on the given data block
                task_params = self.task_function_params
                task_params[self.main_data_param_name] = self.main_data[block_selection]
                self.result.append(self.task_function(**task_params))
                self.blocks.append(block_selection)
                # Record the timings
                end_time = time.time()
                run_time = end_time - start_time
                self.block_times.append(run_time)

        # Return the result
        return self.result, self.blocks


def imports_mpi(python_object):
    """
    Check whether the given class import mpi

    The implementation inspects the source code of the
    analysis to see if MPI is imported by the code.
    """
    import inspect
    import re
    code = inspect.getsource(python_object)
    object_imports_mpi = re.search('import\s+mpi4py', code) is not None
    object_imports_mpi = object_imports_mpi or re.search('from\s+mpi4py\s+import', code) is not None
    object_imports_mpi = object_imports_mpi or re.search('import\s+mpi4py', code) is not None
    object_imports_mpi = object_imports_mpi or re.search('import\s+mpi', code) is not None
    object_imports_mpi = object_imports_mpi or re.search('from\s+mpi\s+import', code) is not None
    object_imports_mpi = object_imports_mpi or re.search('from\s+omsi.shared.mpi_helper\s+import', code) is not None
    object_imports_mpi = object_imports_mpi or re.search('import\s+omsi.shared.mpi_helper', code) is not None
    return object_imports_mpi


def broadcast(data, comm=None, root=0):
    """
    MPI broadcast operation to broadcast data from one rank to all other ranks

    :param data: The data to be gathered
    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    :param root: The rank where the data is sned from
    :return: The data object
    """
    if MPI_AVAILABLE:
        my_comm = comm if comm is not None else MPI.COMM_WORLD
        return my_comm.bcast(data, root)
    else:
        return data


def gather(data, comm=None, root=0):
    """
    MPI gather operation or return a list with just [data,] if MPI is not available

    :param data: The data to be gathered
    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    :param root: The rank where the data should be collected to. Default value is 0
    :return: List of data objects from all the ranks
    """
    if MPI_AVAILABLE:
        my_comm = comm if comm is not None else MPI.COMM_WORLD
        return my_comm.gather(data, root)
    else:
        return [data, ]


def barrier(comm=None):
    """
    MPI barrier operation or no-op when running without MPI

    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    """
    if MPI_AVAILABLE:
        my_comm = comm if comm is not None else MPI.COMM_WORLD
        my_comm.barrier()
    else:
        pass


def get_rank(comm=None):
    """
    Get the current process rank
    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    :return: The integer index of the rank
    """
    rank = 0
    if MPI_AVAILABLE:
        if comm:
            rank = comm.Get_rank()
        else:
            rank = MPI.COMM_WORLD.Get_rank()
    return rank


def get_size(comm=None):
    """
    Get the size of the current communication domain/
    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    :return: The integer index of the rank
    """
    size = 1
    if MPI_AVAILABLE:
        if comm:
            size = comm.Get_size()
        else:
            size = MPI.COMM_WORLD.Get_size()
    return size


def get_comm_world():
    """
    Get MPI.COMM_WORLD
    :return: mpi communicator or None if MPI is not available
    """
    if MPI_AVAILABLE:
        return MPI.COMM_WORLD
    else:
        return None


def is_mpi_available():
    """
    Check if MPI is available. Same as MPI_AVAILABLE
    :return: bool indicating whether MPI is available
    """
    return MPI_AVAILABLE


def mpi_type_from_dtype(dtype):
    """
    Ge the the corresponding MPI type for the given basic numpy dtype

    :param dtype: Basic numpy dtype to be mapped to the MPI type
    :return: The MPI type or None if not found
    """
    if MPI_AVAILABLE:
        try:
            return MPI.__TypeDict__[dtype.char]
        except AttributeError:  # Older versions of mpi4py use _typedict
            return MPI._typedict[dtype.char]
        except KeyError:
            return None
        except:
            return None
    else:
        return None
