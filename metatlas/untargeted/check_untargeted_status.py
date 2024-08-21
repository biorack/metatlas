import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
import argparse

def add_arguments(parser):
    parser.add_argument('--direct_input', type=str, default=None, help='Input project names from command line as a CSV list')
    parser.add_argument('--background_designator', type=str, default="ExCtrl", help='Input control/background sample names from command line as a CSV list')
    parser.add_argument('--print_recent', type=int, default=50, help='Print status for the most recent N projects')

def check_args(args):
    ##### Check if the input arguments are valid
    if args.direct_input:
        args.direct_input = args.direct_input.split(',')
    if args.background_designator:
        args.background_designator = args.background_designator.split(',')
    if args.direct_input is not None:
        args.print_recent = 100000 # Print all projects if direct_input is specified

def main():
    ##### Kick off the script
    start_message = mzm.start_script(script="check_untargeted_status.py")
    print(start_message)

    ##### Set arguments to pass to the pipeline steps and run some checks
    parser = argparse.ArgumentParser(description='Run untargeted pipeline.')
    add_arguments(parser)
    args = parser.parse_args()

    print('Step 1/3: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks(background_designator=args.background_designator)
    print('Step 2/3: Checking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks()
    print('Step 3/3: Getting job status for untargeted projects...')
    mzm.get_untargeted_status(direct_input=args.direct_input,print_recent=args.print_recent)

    ##### End the script
    end_message = mzm.end_script(script="check_untargeted_status.py")
    print(end_message)

if __name__ == "__main__":
    main()