#import sys

from collections import defaultdict
from metatlas.metatlas_objects import find_invalid_runs
from metatlas.system_utils import send_mail

def check_metatlas():
    """
    Checks for invalid runs.
    """
    invalid_runs = find_invalid_runs(_override=True)

    if invalid_runs:
        grouped = defaultdict(list)
        for run in invalid_runs:
            grouped[run.username].append(run.mzml_file)
        for (username, filenames) in grouped.items():
            body = 'You have runs that are not longer accessible\n'
            body += 'To remove them from the database, run the following on ipython.nersc.gov:\n\n'
            body += 'from metatlas.metatlas_objects import find_invalid_runs, remove_objects\n'
            body += 'remove_objects(find_invalid_runs())\n\n'
            body += 'The invalid runs are:\n%s' % ('\n'.join(filenames))
            send_mail('Metatlas Runs are Invalid', username, body)

if __name__ == '__main__':
    check_metatlas()
