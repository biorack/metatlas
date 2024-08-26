import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
import argparse

def add_arguments(parser):
    parser.add_argument('--direct_input', type=str, default=None, help='Input project names from command line as a CSV list')
    parser.add_argument('--background_designator', type=str, default="ExCtrl", help='Input control/background sample names from command line as a CSV list')
    parser.add_argument('--print_recent', type=int, default=50, help='Print status for the most recent N projects')
    parser.add_argument('--update_lims_table', type=bool, default=False, help='Update LIMS table before getting job status')

def check_args(args):
    ##### Check if the input arguments are valid
    if args.direct_input:
        args.direct_input = args.direct_input.split(',')
    if args.background_designator:
        args.background_designator = args.background_designator.split(',')
    if args.direct_input is not None:
        args.print_recent = 100000 # Print all projects if direct_input is specified

def main():

    ##### Set arguments to pass to the pipeline steps and run some checks
    parser = argparse.ArgumentParser(description='Run untargeted pipeline.')
    add_arguments(parser)
    args = parser.parse_args()

    if args.update_lims_table is True:
        print('Checking and updating status of MZmine jobs in LIMS...')
        mzm.update_mzmine_status_in_untargeted_tasks(background_designator=args.background_designator)
        print('Checking and updating status of FBMN jobs in LIMS...')
        mzm.update_fbmn_status_in_untargeted_tasks()
    print('Getting job status for untargeted projects...')
    mzm.get_untargeted_status(direct_input=args.direct_input,print_recent=args.print_recent)

if __name__ == "__main__":
    main()