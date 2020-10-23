# import argparse
from __future__ import absolute_import
from __future__ import print_function
import smtplib
import mimetypes
import itertools
import time
import sys
import os
import glob
import pandas as pd
import metatlas.metatlas_objects.metatlas_objects as metob
import multiprocessing as mp
import numpy as np

from datetime import datetime, time as dtime
from collections import defaultdict
from metatlas.io.mzml_loader import VERSION_TIMESTAMP
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email import encoders
from subprocess import check_output
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.metatlas_objects.metatlas_objects import find_invalid_runs
from metatlas.io.system_utils import send_mail
from six.moves import range
from six.moves import zip


# TO-DO: have these vars be defined from external YML (probably nersc.yml?)

# the lower the tolerance, the more relevant data you might see
# at a cost of being less flexible with defined peaks
# default is 15
tolerance = 15
# the lower the threshold, the more "peaks" you'll see
threshold = 1e5
# change standard if needed. This is ABMBA
std = 229.981116596
std_neg = 227.966564596
std_name = 'ABMBA'
rt_max = -1
run_times = {}
save_path = '/project/projectdirs/metatlas/projects/istd_logs/'


def check_metatlas():
    """
    Checks for invalid runs.
    """
    invalid_runs = find_invalid_runs(_override=True)
    # Currently there exists a bug (?) where it shows empty files.
    # It could be a result of the database remembering invalid files
    # that were removed. Since there is no file associated with these runs,
    # there will not be any filename.
    # This can technically be fixed by running:
    # grouped = grouped.filter(False, str_list)
    # followed by a run to remove all usernames without any files,
    # but does not directly solve the problem.
    # Leaving this to someone in the future to resolve.
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


def yield_ppm(dataset, name, threshold, tolerance, mz_t, rt_max):
    """
    Helper for get_ppm. Given a name, a dataset, intensity threshold,
    ppm tolerance, time of highest mz, and retention time at that point,
    returns a filtered and bounded dataset for ppm under name,
    """
    # bugs out if you divide by mz_t on the left...
    data_filter = (abs(dataset['mz'] - mz_t) * 1e6) <= tolerance * mz_t
    filtered = dataset[data_filter]
    bounded = filtered[filtered.i >= threshold]
    bounded = bounded.loc[:, list(bounded)]
    bounded[name] = (abs(bounded['mz'] - mz_t) * 1e6 / mz_t)
    if not bounded.empty and rt_max != -1:
        data_filter = abs(bounded.rt - rt_max) <= .5
        bounded = bounded[data_filter]
    return bounded


def get_ppms(dataset, tolerance=15, std=229.981116596, std_neg=227.966564596, threshold=1e4,rt_max=rt_max):
    """
    Given some dataset, gets applicable ppm and returns the resulting data.
    Can change various parameters if desired:
    - a different standard than ABMBA
    - tolerances for what is considered to be reasonable
    - a threshold for what is defined as a peak
    """
    # # if neg file, then use std_neg: ms1_pos should be empty if so
    ms1 = 'ms1_pos'
    ms2 = 'ms2_pos'
    mz_t = std
    data = dataset[ms1]
    if data.empty:
        mz_t = std_neg
        ms1 = 'ms1_neg'
        ms2 = 'ms2_neg'        

    # Looking through ms1_pos
    # MS 1
    data = dataset[ms1]

    ms1_max_i = -1

    # get the max i and rt at max i
    filtered = data[(abs(data['mz'] - mz_t) * 1e6) <= tolerance * mz_t]
    tmp = filtered[filtered.i >= threshold]
    tmp = tmp.loc[:, list(tmp)]
    if not tmp.empty:
        max_i_row = tmp['i'].idxmax()
        ms1_max_i = tmp.loc[max_i_row].i
        rt_max = tmp.loc[max_i_row].rt

    d_1 = yield_ppm(data, 'ppm', threshold, tolerance, mz_t, rt_max)

    # Looking through MS2, precursor
    data = dataset[ms2]
    data_filter = (abs(data['precursor_MZ'] - mz_t) * 1e6) <= tolerance * mz_t
    filtered = data[data_filter]
    thresh = filtered[filtered.i >= threshold]
    if rt_max != -1:
        thresh = thresh[abs(thresh.rt - rt_max) <= .5]

    uniques = thresh.drop_duplicates('precursor_MZ')
    d_2 = uniques.loc[:, list(uniques)]
    d_2['ppm'] = (abs(d_2['precursor_MZ'] - mz_t) * 1e6 / mz_t)

    # Double filter on mz and precursor_MZ
    # MS 2
    d_3 = yield_ppm(filtered, 'ppm', threshold, tolerance, mz_t, rt_max)

    return (d_1.reset_index(), d_2.reset_index(), d_3.reset_index(),
            ms1_max_i, rt_max)


def data_verify(file_name,tolerance=tolerance,std=std,std_neg=std_neg,threshold=threshold,rt_max=-1):
    if type(file_name) == str:
        hdf5_name = file_name
    else:
        #its a metatlas object
        hdf5_name = file_name.hdf5_file

    dataset = ma_data.df_container_from_metatlas_file(hdf5_name)
    s = hdf5_name.split('/')
    samples_dict = {}
    ppms = get_ppms(dataset, tolerance=tolerance, std=std,std_neg=std_neg, threshold=threshold,rt_max=rt_max)
    samples_dict['file name'] = s[-1]
    samples_dict['ms1 ppm'] = ppms[0]['ppm'].tolist()
    samples_dict['ms1 intensity'] = ppms[0]['i'].tolist()
    samples_dict['ms2_precursor ppm'] = ppms[1]['ppm'].tolist()
    samples_dict['ms2 ppm'] = ppms[2]['ppm'].tolist()

    samples_dict['rt'] = ppms[0]['rt'].tolist()

    # extra stats and more data filtering
    samples_dict['ms1 ppm weighted mean'] = weighted_mean(ppms[0])
    # hacky
    samples_dict['ms2_precursor ppm weighted mean'] = weighted_mean(ppms[1])
    samples_dict['ms2 ppm weighted mean'] = weighted_mean(ppms[2])

    samples_dict['ms1 max intensity'] = ppms[3]
    samples_dict['rt at max intensity'] = ppms[4]
    if type(file_name) == str:
        samples_dict['acquisition timestamp'] = run_times[u'/global%s' % file_name]
    else:
        #its a metatlas object
        samples_dict['acquisition timestamp'] = file_name.acquisition_time
        
    samples_dict['acquisition time'] = datetime.fromtimestamp(
        samples_dict['acquisition timestamp']).strftime('%Y-%m-%d %H:%M:%S')

    for key in samples_dict.keys():
        if samples_dict[key] == [] or samples_dict[key] == -1:
            samples_dict[key] = 0

    return samples_dict


def process_files(files):
    """
    Goes through the files and returns a list of dicts which contain info about
    the runs. Returns -1 if it failed.
    """
    samples = []
    if files:
        try:
            print('found %d files' % len(files))
            p = mp.Pool(min(32, len(files)))
            samples = p.map(data_verify, files)
        finally:
            p.close()
            p.join()

            print('done processing files')
        return samples
    return -1


def clean_up(samples, files):
    """
    Clean up on a bunch of files to generate CSVs and a report.
    """
    directory = files[0].split('/')
    path = directory[len(directory) - 2]
    df = pd.DataFrame(samples)
    sorted_cols = ['file name', 'ms1 ppm', 'ms1 ppm weighted mean',
                   'ms1 intensity', 'ms1 max intensity']
    sorted_cols.extend(['rt', 'rt at max intensity', 'ms2_precursor ppm',
                        'ms2_precursor ppm weighted mean'])
    sorted_cols.extend(['ms2 ppm', 'ms2 ppm weighted mean',
                        'acquisition timestamp', 'acquisition time'])

    # get the full report
    df = df.loc[:, sorted_cols]
    df = df.sort_values('file name')
    csv_name = save_path + '%s_run.csv' % path
    df.to_csv(csv_name, index=False)

    print('begin clean up')
    # abridged report
    df = df.sort_values('file name')
    sorted_cols = ['file name', 'acquisition time', 'ms1 ppm weighted mean',
                   'ms1 max intensity', 'rt at max intensity']
    sorted_cols.extend(['ms2_precursor ppm weighted mean',
                        'ms2 ppm weighted mean'])
    df = df.loc[:, sorted_cols]
    for col in sorted_cols:
        if 'weighted mean' in col or 'rt at max' in col:
            for i in range(df.shape[0]):
                if df.loc[i, col] != 0:
                    df.loc[i, col] = '%.2f' % df[col][i]
        elif col == 'ms1 max intensity':
            for i in range(df.shape[0]):
                if df.loc[i, col] != 0:
                    df.loc[i, col] = '%.3e' % df[col][i]
    data_filter = df['file name'].str.contains('Blank')
    df_blanks = df[data_filter]
    df_non_blanks = df[~data_filter]
    print('done clean up')
    return [csv_name, (df_non_blanks, df_blanks)]


def weighted_mean(df):
    """
    Calculates the weighted mean for some ppm and intensity data.
    """
    i_sum = df['i'].sum()
    if i_sum:
        vm = df['i'].dot(df['ppm'])
        return vm / i_sum
    return -1


def convert_to_html(df):
    """
    Converts a dataframe into a styled html table for emails.
    """
    pd.options.display.max_colwidth = -1
    html = df.to_html(index=False)
    # hacky way to get css in the file
    html = html.replace('<table border="1"',
                        '<table style="width: 100%;')
    html = html.replace('<th>',
                        '<th style="border-bottom: 1px solid #ddd; padding: 8px; text-align: left; margin: 0">')
    html = html.replace('<td>',
                        '<td style="border-bottom: 1px solid #ddd; padding: 8px; text-align: left; margin: 0">')
    html = html.replace('<tr style="text-align: right;">',
                        '<tr style="text-align: left">')
    return html


def send_run_email(username, run_file, tables, standard=std_name):
    """
    Sends an email of the run(s) done and packages the generated CSV(s).
    """
    sender = 'pasteur@nersc.gov'
    receivers = ['%s@nersc.gov' % username, '%s@nersc.gov' % 'bpb']
    # Build email
    msg = MIMEMultipart()
    msg['From'] = sender
    msg['To'] = ", ".join(receivers)  # needs to be a string
    msg['Subject'] = "Run Information"

    template = "Hello %s,\n\n" % username
    template += "Here is a daily summary of your runs with %s from the last 3 days.\n" % standard
    template += "0's in the table indicate that the particular value was not detected.\n"
    template += "Attached is a more verbose summary of the run in a .csv file.\n"
    # add more to the template if needed?

    text = MIMEText(template, 'plain')
    msg.attach(text)

    for title, table in zip(run_file, tables):
        s = title.replace('run_data/', '').replace('_run.csv', '\n')
        num_files = len(table[0]) + len(table[1])
        txt = '<h2 style="border-top: solid; padding-top: 10px;">Run Report for %s:</h2>' % s
        txt += '<p>Detected %d runs!</p>' % num_files
        text = MIMEText(txt, 'html')
        msg.attach(text)
        table_non_blank = convert_to_html(table[0])
        html = MIMEText(table_non_blank, 'html')
        msg.attach(html)
        text = MIMEText('<h3>Blanks in %s:</h3>' % s, 'html')
        msg.attach(text)
        table_blanks = convert_to_html(table[1])
        html = MIMEText(convert_to_html(table[1]), 'html')
        msg.attach(html)

    for f in run_file:
        ctype, encoding = mimetypes.guess_type(f)
        if ctype is None or encoding is not None:
            ctype = "application/octet-stream"

        maintype, subtype = ctype.split("/", 1)
        # Parser for generic files found on stackoverflow
        maintype, subtype = ctype.split("/", 1)
        fp = open(f, "rb")
        attachment = MIMEBase(maintype, subtype)
        attachment.set_payload(fp.read())
        fp.close()
        encoders.encode_base64(attachment)
        # another case where save path needs to be modular / in a .yml
        attachment.add_header("Content-Disposition", "attachment",
                              filename=f.replace(save_path, ''))
        msg.attach(attachment)

    try:
        smtpObj = smtplib.SMTP('localhost')
        smtpObj.sendmail(sender, receivers, msg.as_string())
        sys.stdout.write("Successfully sent email to %s\n" % username)
        sys.stdout.flush()
    except smtplib.SMTPException:
        sys.stderr.write("Error: unable to send email to %s\n" % username)
        sys.stdout.flush()


# Run this to start the task
def run_checker():
    # global std
    # global std_neg
    # global std_name
    # parser = argparse.ArgumentParser(description="Parameters for custom compounds")
    # parser.add_argument('-mz_pos', '--mz_pos', type=float, required=False)
    # parser.add_argument('-mz_neg', '--mz_neg', type=float, required=False)
    
    # args = parser.parse_args()
    # if args.mz_pos is not None:
    #     std = args.mz_pos
    #     std_neg = args.mz_neg

    print('starting task')
    start_time = time.time()
    origin_directory = '/project/projectdirs/metatlas/raw_data/'
    # pull users, email them their run info
    delta = 60*24*3  # 3 days to minutes
    check = 'find %s -name "*.h5" -mmin -%d -size +2k' % (
            origin_directory, delta)
    valid_files = check_output(check, shell=True).decode('utf-8').splitlines()
    valid_files = set(valid_files)

    # dict of user name and the directory of their files
    files = defaultdict(set)
    for i in valid_files:
        s = i.split('/')
        usr = s[-3]
        directory = s[-2]
        files[usr].add(directory)

    print('detected a total of %d file(s)!\n' % len(valid_files))

    for username in files.keys():
        print('observing files from', username)
        csvs = []
        html_tables = []
        for directory in files[username]:
            print('looking at', directory)
            extension = '%s/%s/' % (username, directory)
            # list of all files in a given directory
            dir_files = glob.glob(os.path.join(origin_directory +
                                  extension, '*.h5'))
            # contains a list of dicts
            samples = process_files(dir_files)
            # a list containing a csv name
            # and a tuple of the non-blank, blank dataframes
            info = clean_up(samples, dir_files)
            csvs.append(info[0])
            html_tables.append((info[1][0], info[1][1]))
        print('attempting to send email')
        send_run_email(username, csvs, html_tables)
        print('\n')  # for formatting purposes

    print('finished task in %s seconds' % (time.time() - start_time))

# Integrity monitor does two major things:
# - Checks the database for invalid files and warns the user.
# - Runs through all the raw data collected in the past 3 days and gives the
#   user a general report of the data collected.
if __name__ == '__main__':
    check_metatlas()
    runs = metob.retrieve('lcmsruns', username='*')
    for i in runs:
        run_times[i.hdf5_file] = i.acquisition_time
    run_checker()
