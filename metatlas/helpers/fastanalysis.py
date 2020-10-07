import sys
import os
import multiprocessing as mp
import pprint

from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
from metatlas import metatlas_objects as metob
from metatlas.helpers import dill2plots as dp
from metatlas.helpers import chromplotplus as cpp
from metatlas.helpers import spectralprocessing as sp

import numpy as np
import pandas as pd

loose_param = {'min_intensity': 1e3,
               'rt_tolerance': .25,
               'mz_tolerance': 25,
               'min_msms_score': 0.3, 'allow_no_msms': True,
               'min_num_frag_matches': 1,  'min_relative_frag_intensity': .01}

strict_param = {'min_intensity': 1e5,
                'rt_tolerance': .25,
                'mz_tolerance': 5,
                'min_msms_score': .6, 'allow_no_msms': False,
                'min_num_frag_matches': 3, 'min_relative_frag_intensity': .1}

def make_stats_table(input_fname = '', input_dataset = [], msms_hits_df = None,
                     include_lcmsruns = [], exclude_lcmsruns = [], include_groups = [], exclude_groups = [],
                     output_loc = None,
                     output_sheetname = 'Draft_Final_Idenfications.xlsx',
                     msms_hits = None,
                     min_peak_height=0, min_num_data_points=0, rt_tolerance=np.inf, ppm_tolerance=np.inf,

                     min_msms_score=0, min_num_frag_matches=0,
                     allow_no_msms=False, min_relative_frag_intensity=None,
                     use_labels=False, return_all=False,
                     msms_refs_loc='/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v2.tab',
                     dependencies = {'peak_height': [],
                                     'peak_area': ['peak_height'],
                                     'num_data_points': ['peak_height'],
                                     'rt_peak': ['peak_height', 'rt_delta'],
                                     'rt_delta': ['peak_height'],
                                     'mz_centroid': ['peak_height', 'mz_ppm'],
                                     'mz_ppm': ['peak_height'],
                                     'msms_score': ['peak_height', 'num_frag_matches'],
                                     'num_frag_matches': ['peak_height', 'msms_score']}):

    assert output_loc is not None or return_all

    if not input_dataset:
        metatlas_dataset = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        metatlas_dataset = input_dataset

    if output_loc is not None and not os.path.exists(output_loc):
        os.mkdir(output_loc)
    if output_loc is not None and not os.path.exists(os.path.join(output_loc,'data_sheets')):
        os.mkdir(os.path.join(output_loc,'data_sheets'))
    if output_loc is not None and not os.path.exists(os.path.join(output_loc,'stats_tables')):
        os.mkdir(os.path.join(output_loc,'stats_tables'))

    # filter runs from the metatlas dataset
    if include_lcmsruns:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_include_list(metatlas_dataset,'lcmsrun',include_lcmsruns)
    if include_groups:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_include_list(metatlas_dataset,'group',include_groups)

    if exclude_lcmsruns:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_exclude_list(metatlas_dataset,'lcmsrun',exclude_lcmsruns)
    if exclude_groups:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_exclude_list(metatlas_dataset,'group',exclude_groups)

    final_df = pd.DataFrame(columns=['index'])
    file_names = ma_data.get_file_names(metatlas_dataset)
    compound_names = ma_data.get_compound_names(metatlas_dataset,use_labels=use_labels)[0]

    metrics = ['msms_score', 'num_frag_matches', 'mz_centroid', 'mz_ppm', 'rt_peak', 'rt_delta', 'peak_height', 'peak_area', 'num_data_points']

    dfs = {m:None for m in metrics}
    passing = {m:np.ones((len(compound_names), len(file_names))).astype(float) for m in metrics}

    for metric in ['peak_height', 'peak_area', 'rt_peak', 'mz_centroid']:
        dfs[metric] = dp.make_output_dataframe(input_dataset=metatlas_dataset, fieldname=metric, use_labels=use_labels)

    dfs['mz_ppm'] = dfs['peak_height'].copy()
    dfs['mz_ppm'] *= np.nan

    dfs['num_data_points'] = pd.DataFrame([[len(metatlas_dataset[i][j]['data']['eic']['intensity'])
                                           for i in range(len(metatlas_dataset))]
                                          for j in range(len(metatlas_dataset[0]))])
    dfs['num_data_points'].index = dfs['mz_ppm'].index
    dfs['msms_score'] = dfs['mz_ppm'].copy()
    dfs['num_frag_matches'] = dfs['mz_ppm'].copy()
    dfs['rt_delta'] = dfs['mz_ppm'].copy()

    passing['peak_height'] = (np.nan_to_num(dfs['peak_height'].values) >= min_peak_height).astype(float)
    passing['num_data_points'] = (np.nan_to_num(dfs['num_data_points'].values) >= min_num_data_points).astype(float)

    #msms_hits_df = dp.get_msms_hits(metatlas_dataset, use_labels, ref_index=['database', 'id', 'inchi_key', 'precursor_mz'])
    #msms_hits_df = dp.get_msms_hits(metatlas_dataset, use_labels, ref_index=['database', 'id', 'inchi_key'])
    #msms_hits_df.rename(columns={'inchi_key':'inchi_key_2'},inplace=True)
    msms_hits_df = msms_hits.copy()
    msms_hits_df.reset_index(inplace=True)

    for compound_idx, compound_name in enumerate(compound_names):

        ref_rt_peak = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_peak
        ref_mz = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz

        dfs['rt_delta'].iloc[compound_idx] = abs(ref_rt_peak - dfs['rt_peak'].iloc[compound_idx])
        passing['rt_delta'][compound_idx] = (abs(ref_rt_peak - np.nan_to_num(dfs['rt_peak'].iloc[compound_idx].values)) <= rt_tolerance).astype(float)

        dfs['mz_ppm'].iloc[compound_idx] = 1e6*(abs(ref_mz - dfs['mz_centroid'].iloc[compound_idx]) / ref_mz)
        passing['mz_ppm'][compound_idx] = (dfs['mz_ppm'].iloc[compound_idx].values <= ppm_tolerance).astype(float)

        try:
            inchi_key = metatlas_dataset[0][compound_idx]['identification'].compound[0].inchi_key
        except:
            inchi_key = ''
        compound_ref_rt_min = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_min
        compound_ref_rt_max = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_max
        cid = metatlas_dataset[0][compound_idx]['identification']
        mz_theoretical = cid.mz_references[0].mz
        mz_measured = metatlas_dataset[0][compound_idx]['data']['ms1_summary']['mz_centroid']
        delta_mz = abs(mz_theoretical - mz_measured)
        delta_ppm = delta_mz / mz_theoretical * 1e6

        comp_msms_hits = msms_hits_df[(msms_hits_df['inchi_key'] == inchi_key) \
                                    & (msms_hits_df['msms_scan'] >= compound_ref_rt_min) \
                                    & (msms_hits_df['msms_scan'] <= compound_ref_rt_max) \
                                    & ((abs(msms_hits_df['measured_precursor_mz'].values.astype(float) - mz_theoretical)/mz_theoretical) \
                                    <= cid.mz_references[0].mz_tolerance*1e-6)]

        comp_msms_hits = comp_msms_hits.sort_values('score', ascending=False)
        file_idxs, scores, msv_sample_list, msv_ref_list, rt_list = [], [], [], [], []
        if len(comp_msms_hits) > 0 and not np.isnan(np.concatenate(comp_msms_hits['msv_ref_aligned'].values, axis=1)).all():
            file_idxs = [file_names.index(f) for f in comp_msms_hits['file_name']]
            scores = comp_msms_hits['score'].values.tolist()
            msv_sample_list = comp_msms_hits['msv_query_aligned'].values.tolist()
            msv_ref_list = comp_msms_hits['msv_ref_aligned'].values.tolist()
            rt_list = comp_msms_hits['msms_scan'].values.tolist()
            mz_sample_matches = sp.partition_aligned_ms_vectors(msv_sample_list[0], msv_ref_list[0])[0][0].tolist()

        avg_mz_measured = []
        avg_rt_measured = []
        intensities = pd.DataFrame()
        for file_idx, file_name in enumerate(file_names):
            if not np.isnan(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['mz_centroid']):
                avg_mz_measured.append(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['mz_centroid'])
            if not np.isnan(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['rt_peak']):
                avg_rt_measured.append(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['rt_peak'])
            if not np.isnan(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['peak_height']):
                intensities.loc[file_idx, 'file_id'] = file_idx
                intensities.loc[file_idx, 'intensity'] = metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['peak_height']

        avg_mz_measured = np.mean(avg_mz_measured)
        avg_rt_measured = np.mean(avg_rt_measured)

        final_df = final_df.append({'index':compound_idx}, ignore_index=True)
        final_df.loc[compound_idx, 'identified_metabolite'] = ""
        if use_labels or len(cid.compound) == 0:
            cid_label = cid.name
            final_df.loc[compound_idx, 'label'] = cid_label
        else:
            cid_label = cid.compound[0].name
            final_df.loc[compound_idx, 'label'] = cid_label
        
        overlapping_compounds = []
        inchi_key_map = {}

        if(len(cid.compound) != 0):
            #Loop through compounds to identify overlapping compounds
            for compound_iterator in range(len(compound_names)):
                if len(metatlas_dataset[0][compound_iterator]['identification'].compound) == 0:
                    continue
                if use_labels:
                    cpd_iter_label = metatlas_dataset[0][compound_iterator]['identification'].name
                else:
                    cpd_iter_label = metatlas_dataset[0][compound_iterator]['identification'].compound[0].name
                cpd_iter_id = metatlas_dataset[0][compound_iterator]['identification']
                cpd_iter_mz = cpd_iter_id.mz_references[0].mz
                cid_mass = cid.compound[0].mono_isotopic_molecular_weight
                cpd_iter_mass = cpd_iter_id.compound[0].mono_isotopic_molecular_weight
                cid_rt_min = cid.rt_references[0].rt_min
                cid_rt_max = cid.rt_references[0].rt_max
                cpd_iter_rt_min = cpd_iter_id.rt_references[0].rt_min
                cpd_iter_rt_max = cpd_iter_id.rt_references[0].rt_max
                if compound_idx != compound_iterator:
                    if ((cpd_iter_mz-0.005 <= mz_theoretical <= cpd_iter_mz+0.005) or (cpd_iter_mass-0.005 <= cid_mass <= cpd_iter_mass+0.005)) and \
                            ((cpd_iter_rt_min <= cid_rt_min <=cpd_iter_rt_max) or (cpd_iter_rt_min <= cid_rt_max <= cpd_iter_rt_max) or \
                            (cid_rt_min <= cpd_iter_rt_min <= cid_rt_max) or (cid_rt_min <= cpd_iter_rt_max <= cid_rt_max)):
                        overlapping_compounds.append(cpd_iter_label)
                        inchi_key_map[cpd_iter_label] = cpd_iter_id.compound[0].inchi_key

        if len(overlapping_compounds) > 0:
            overlapping_compounds.append(cid_label)
            if len(cid.compound) > 0:
                inchi_key_map[cid_label] = cid.compound[0].inchi_key
            else:
                inchi_key_map[cid_label] = ""
            final_df.loc[compound_idx, 'overlapping_compound'] = "//".join(cpd for cpd in sorted(overlapping_compounds, key=str))
            final_df.loc[compound_idx, 'overlapping_inchi_keys'] = "//".join(inchi_key_map[cpd] for cpd in sorted(overlapping_compounds, key=str))
        else:
            final_df.loc[compound_idx, 'overlapping_compound'] = ""
            final_df.loc[compound_idx, 'overlapping_inchi_keys'] = ""
        if len(cid.compound) == 0:
            final_df.loc[compound_idx, 'formula'] = ""
            final_df.loc[compound_idx, 'polarity'] = cid.mz_references[0].detected_polarity
            final_df.loc[compound_idx, 'exact_mass'] = ""
            final_df.loc[compound_idx, 'inchi_key'] = ""
        else:
            final_df.loc[compound_idx, 'formula'] = cid.compound[0].formula
            final_df.loc[compound_idx, 'polarity'] = cid.mz_references[0].detected_polarity
            final_df.loc[compound_idx, 'exact_mass'] = cid.compound[0].mono_isotopic_molecular_weight
            final_df.loc[compound_idx, 'inchi_key'] = cid.compound[0].inchi_key
        final_df.loc[compound_idx, 'msms_quality'] = ""
        final_df.loc[compound_idx, 'mz_quality'] = ""
        final_df.loc[compound_idx, 'rt_quality'] = ""
        final_df.loc[compound_idx, 'total_score'] = ""
        final_df.loc[compound_idx, 'msi_level'] = ""
        final_df.loc[compound_idx, 'isomer_details'] = ""
        final_df.loc[compound_idx, 'identification_notes'] = cid.identification_notes
        final_df.loc[compound_idx, 'ms1_notes'] = cid.ms1_notes
        final_df.loc[compound_idx, 'ms2_notes'] = cid.ms2_notes
        if len(intensities) > 0:
            final_df.loc[compound_idx, 'max_intensity'] = intensities.loc[intensities['intensity'].idxmax()]['intensity']
            max_intensity_file_id = int(intensities.loc[intensities['intensity'].idxmax()]['file_id'])
            final_df.loc[compound_idx, 'max_intensity_file'] = file_names[max_intensity_file_id]
            final_df.loc[compound_idx, 'ms1_rt_peak'] = metatlas_dataset[max_intensity_file_id][compound_idx]['identification'].rt_references[0].rt_peak
        else:
            final_df.loc[compound_idx, 'max_intensity'] = ""
            final_df.loc[compound_idx, 'max_intensity_file'] = ""
            final_df.loc[compound_idx, 'ms1_rt_peak'] = ""
        if file_idxs != []:
            final_df.loc[compound_idx, 'msms_file'] = file_names[file_idxs[0]]
            final_df.loc[compound_idx, 'msms_rt'] = float("%.2f" % rt_list[0])
            final_df.loc[compound_idx, 'msms_numberofions'] = len(mz_sample_matches)
            final_df.loc[compound_idx, 'msms_matchingions'] = ','.join(['%5.3f'%m for m in mz_sample_matches])
            if len(mz_sample_matches) == 1:
                # Set score to zero when there is only one matching ion. precursor intensity is set as score in such cases and need to be set to 0 for final identification.
                final_df.loc[compound_idx, 'msms_score'] = 0.0
            else:
                final_df.loc[compound_idx, 'msms_score'] = float("%.4f" % scores[0])
        else:
            final_df.loc[compound_idx, 'msms_file'] = ""
            final_df.loc[compound_idx, 'msms_rt'] = ""
            final_df.loc[compound_idx, 'msms_numberofions'] = ""
            final_df.loc[compound_idx, 'msms_matchingions'] = ""
            final_df.loc[compound_idx, 'msms_score'] = ""
        final_df.loc[compound_idx, 'mz_adduct'] = cid.mz_references[0].adduct
        final_df.loc[compound_idx, 'mz_theoretical'] = float("%.4f" % mz_theoretical)
        final_df.loc[compound_idx, 'mz_measured'] = float("%.4f" % avg_mz_measured)
        final_df.loc[compound_idx, 'mz_error'] = float("%.4f" % abs(mz_theoretical - avg_mz_measured))
        final_df.loc[compound_idx, 'mz_ppmerror'] = float("%.4f" % (abs(mz_theoretical - avg_mz_measured) / mz_theoretical * 1e6))
        final_df.loc[compound_idx, 'rt_min'] = float("%.2f" % compound_ref_rt_min)
        final_df.loc[compound_idx, 'rt_max'] = float("%.2f" % compound_ref_rt_max)
        final_df.loc[compound_idx, 'rt_theoretical'] = float("%.2f" % cid.rt_references[0].rt_peak)
        final_df.loc[compound_idx, 'rt_measured'] = float("%.2f" % avg_rt_measured)
        final_df.loc[compound_idx, 'rt_error'] = float("%.2f" % abs(cid.rt_references[0].rt_peak - avg_rt_measured))


        for file_idx, file_name in enumerate(file_names):
            if len(msms_hits_df) == 0:
                rows = []
            else:
                rows = msms_hits_df[(msms_hits_df['inchi_key'] == inchi_key) & \
                                (msms_hits_df['file_name'] == file_name) & \
                                (msms_hits_df['msms_scan'] >= compound_ref_rt_min) & (msms_hits_df['msms_scan'] <= compound_ref_rt_max) & \
                                ((abs(msms_hits_df['measured_precursor_mz'].values.astype(float) - metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz)/metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz) \
                                   <= metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz_tolerance*1e-6)]

            if len(rows) == 0:
                dfs['msms_score'].iat[compound_idx, file_idx] = np.nan
                dfs['num_frag_matches'].iat[compound_idx, file_idx] = np.nan
            else:
                if not np.isnan(np.concatenate(rows['msv_ref_aligned'].values, axis=1)).all():
                    dfs['msms_score'].iat[compound_idx, file_idx] = rows.loc[rows['score'].idxmax()]['score']
                dfs['num_frag_matches'].iat[compound_idx, file_idx] = rows.loc[rows['score'].idxmax()]['num_matches']

    passing['msms_score'] = (np.nan_to_num(dfs['msms_score'].values) >= min_msms_score).astype(float)
    passing['num_frag_matches'] = (np.nan_to_num(dfs['num_frag_matches'].values) >= min_num_frag_matches).astype(float)

    if not output_sheetname.endswith('.xlsx'):
        output_sheetname = output_sheetname + '.xlsx'
    writer = pd.ExcelWriter(os.path.join(output_loc,output_sheetname), engine='xlsxwriter')
    final_df.to_excel(writer, sheet_name='Final_Identifications', index=False, startrow=3)
    
    #set format
    workbook = writer.book
    f_blue = workbook.add_format({'bg_color': '#DCFFFF'})
    f_yellow = workbook.add_format({'bg_color': '#FFFFDC'})
    f_rose = workbook.add_format({'bg_color': '#FFDCFF'})
    cell_format = workbook.add_format({'bold': True, 'align': 'center'})
    scientific_format = workbook.add_format({'num_format': '0.00E+00'})
    cell_format.set_text_wrap()
    cell_format.set_border()
    worksheet = writer.sheets['Final_Identifications']
    worksheet.set_row(1,60)
    worksheet.set_row(2,60)
    worksheet.set_column('S:S', None, scientific_format)
    worksheet.merge_range('A1:I1', 'COMPOUND ANNOTATION', cell_format)
    worksheet.merge_range('J1:R1', 'COMPOUND IDENTIFICATION SCORES', cell_format)
    worksheet.merge_range('S1:U1', 'MS1 INTENSITY INFORMATION', cell_format)
    worksheet.merge_range('V1:Y1', 'MSMS INFORMATION', cell_format)
    worksheet.write('Z1', 'MSMS EVALUATION', cell_format)
    worksheet.merge_range('AA1:AC1', 'ION INFORMATION', cell_format)
    worksheet.merge_range('AD1:AE1', 'M/Z EVALUATION', cell_format)
    worksheet.merge_range('AF1:AI1', 'CHROMATOGRAPHIC PEAK INFORMATION', cell_format)
    worksheet.write('AJ1', 'RT EVALUATION', cell_format)

    HEADER2 = ['Compound #','Identified Metabolite','Name of metabolite searched for','Labels of Overlapping Compounds','Inchi Keys of Overlapping Compounds','Molecular Formula','Polarity','Exact Mass','Inchi Key','MSMS Score (0 to 1)','m/z score (0 to 1)','RT score (0 to 1)','Total ID Score (0 to 3)','Mass Spec Inititative Identification Level','Isomer details','Identification notes','MS1 notes', 'MS2 notes','Maximum MS1 intensity across all files','Filename w/ maximum MS1','Retention time of max intensity MS1 peak','File with highest MSMS match score','RT of highest matched MSMS scan','Number of ion matches in msms spectra to EMA reference spectra','List of ion matches in msms spectra to EMA reference spectra','','Adduct','Theoretical m/z','Measured m/z','mass error (delta Da)','mass error (delta ppm)','Minimum retention time (min)','Maximum retention time (max)','Theoretical retention time (peak)','Detected RT (peak)','RT error (absolute delta)']

    for i, header in enumerate(HEADER2):
        worksheet.write(1,i, header, cell_format)

    HEADER3 = ['Unique for study','Some isomers are not chromatographically or spectrally resolvable.','Name of standard reference compound in library match.','compound with similar mz (abs difference <= 0.005) or monoisotopic molecular weight (abs difference <= 0.005) and RT (min or max within the RT-min-max-range of similar compound)','List of inchi keys that correspond to the compounds listed in the previous column','','','monoisotopic mass (neutral except for permanently charged molecules)','neutralized version','1 (MSMS matches ref. std.), 0.5 (possible match), 0 (no MSMS collected or no appropriate ref available), -1 (bad match)','POS MODE: 1 (delta ppm </= 5 or delta Da </= 0.001), 0.5 (delta ppm 5-10 and delta Da > 0.001), 0 (delta ppm > 10) NEG MODE: 1 (delta ppm </= 15 or delta Da </= 0.005), 0.5 (delta ppm 15-20 and delta Da > 0.001), 0 (delta ppm > 20) (NOTE: neg mass accuracy is not as good as pos) mz_quality','1 (delta RT </= 0.5), 0.5 (delta RT > 0.5 & </= 2), 0 (delta RT > 2 min)','sum of m/z, RT and MSMS score','Level 1 = Two independent and orthogonal properties match authentic standard; else = putative [Metabolomics. 2007 Sep; 3(3): 211-221. doi: 10.1007/s11306-007-0082-2]','Isomers have same formula (and m/z) and similar RT - MSMS spectra may be used to differentiate (exceptions) or RT elution order','','','','','','','','','mean # of fragment ions matching between compound in sample and reference compound / standard; may include parent and isotope ions and very low intensity background ions (these do not contribute to score)','','MSMS score (highest across all samples), scale of 0 to 1 based on an algorithm. 0 = no match, 1 = perfect match. If no score, then no MSMS was acquired for that compound (@ m/z & RT window).','More than one may be detectable; the one evaluated is listed','theoretical m/z for a given compound / adduct pair','average m/z within 20ppm of theoretical detected across all samples @ RT peak','absolute difference between theoretical and detected m/z','ppm difference between theoretical and detected m/z','','','theoretical retention time for a compound based upon reference standard at highest intensity point of peak','average retention time for a detected compound at highest intensity point of peak across all samples','absolute difference between theoretical and detected RT peak']
    
    for i, header in enumerate(HEADER3):
        worksheet.write(2,i, header, cell_format)

    worksheet.merge_range('AF3:AG3', 'Retention range including start and end of detection of an m/z value (Note: Peak Height is calculated as the highest intensity of an m/z within the min/max RT range. Peak Area is calculated as the integrated area under the curve for an m/z within the mix/max RT range.)', cell_format)
    worksheet.conditional_format('J1:R'+str(len(final_df)+4),{ 'type':'no_errors', 'format':f_blue})
    worksheet.conditional_format('S1:Y'+str(len(final_df)+4),{ 'type':'no_errors', 'format':f_yellow})
    worksheet.conditional_format('Z1:Z'+str(len(final_df)+4),{ 'type':'no_errors', 'format':f_rose})
    worksheet.conditional_format('AA1:AC'+str(len(final_df)+4),{ 'type':'no_errors', 'format':f_yellow})
    worksheet.conditional_format('AD1:AE'+str(len(final_df)+4),{ 'type':'no_errors', 'format':f_rose})
    worksheet.conditional_format('AF1:AI'+str(len(final_df)+4),{ 'type':'no_errors', 'format':f_yellow})
    worksheet.conditional_format('AJ1:AJ'+str(len(final_df)+4),{ 'type':'no_errors', 'format':f_rose})
    writer.save()
    

    #final_df.to_csv(os.path.join(output_loc, 'Draft_Final_Idenfications.tab'), sep='\t')
    for metric in metrics:
        passing[metric][passing[metric] == 0] = np.nan
    stats_table = []


#    for metric in metrics:
#        test = np.product(np.array([passing[dep] for dep in dependencies[metric]]), axis=0)
#        # group_df = (dfs[metric] * test).T.groupby('group').describe()
#        #if output_loc is not None:
#        #    (dfs[metric] * test).to_csv(os.path.join(output_loc, 'data_sheets', 'filtered_%s.tab'%metric), sep='\t')
#        stats_df = (dfs[metric] * test * passing[metric]).T.describe().T
#        stats_df['range'] = stats_df['max'] - stats_df['min']
#        stats_df.columns = pd.MultiIndex.from_product([['filtered'], [metric], stats_df.columns])
#        stats_table.append(stats_df)

    for metric in metrics:
        #if output_loc is not None:
        #    dfs[metric].to_csv(os.path.join(output_loc, 'data_sheets', 'unfiltered_%s.tab'%metric), sep='\t')
        stats_df = dfs[metric].T.describe().T
        stats_df['range'] = stats_df['max'] - stats_df['min']
        #stats_df.columns = pd.MultiIndex.from_product([['unfiltered'], [metric], stats_df.columns])
        stats_df.columns = pd.MultiIndex.from_product([[metric], stats_df.columns])
        stats_table.append(stats_df)

    stats_table = pd.concat(stats_table, axis=1)

    if output_loc is not None:
        stats_table.to_csv(os.path.join(output_loc, 'stats_tables', 'stats_table.tab'), sep='\t')

        with open(os.path.join(output_loc, 'stats_tables/', 'stats_table.readme'), 'w') as readme:
            for var in ['dependencies', 'min_peak_height', 'rt_tolerance', 'ppm_tolerance', 'min_msms_score', 'min_num_frag_matches']:
                readme.write('%s\n'%var)
                try:
                    if np.isinf(eval(var)):
                        pprint.pprint('default', readme)
                    else:
                        pprint.pprint(eval(var), readme)
                except TypeError:
                    pprint.pprint(eval(var), readme)
                readme.write('\n')

    if return_all:
        return stats_table, dfs, passing


def make_scores_df(metatlas_dataset, msms_hits):
    """
    Returns pandas dataframe with columns 'max_intensity', 'median_rt_shift','median_mz_ppm', 'max_msms_score',
    'num_frag_matches', and 'max_relative_frag_intensity', rows of compounds in metatlas_dataset, and values
    of the best "score" for a given compound across all files.

    'max_intensity': highest intensity across all files for given compound
    'median_rt_shift': median shift of RT across all files for given compound to reference
    'median_mz_ppm': median ppm of mz across all files for given compound relative to reference
    'max_msms_score': highest compound dot-product score across all files for given compound relative to reference
    'num_frag_matches': number of matching mzs when calculating max_msms_score
    'max_relative_frag_intensity': ratio of second highest to first highest intensity of matching sample mzs

    :param metatlas_dataset:
    :param msms_hits:

    :return scores_df: pandas dataframe
    """

    file_names = ma_data.get_file_names(metatlas_dataset)
    compound_names = ma_data.get_compound_names(metatlas_dataset)[0]

    scores = []

    #msms_hits_df = dp.get_msms_hits(metatlas_dataset, ref_index=['database', 'id', 'inchi_key', 'precursor_mz'])
    #msms_hits_df = dp.get_msms_hits(metatlas_dataset, ref_index=['database', 'id', 'inchi_key'])
    msms_hits_df = msms_hits.copy()
    #msms_hits_df.rename(columns={'inchi_key':'inchi_key_2'},inplace=True)
    msms_hits_df.reset_index(inplace=True)

    for compound_idx in range(len(compound_names)):
        intensities = []
        rt_shifts = []
        mz_ppms = []
        max_msms_score = np.nan
        num_frag_matches = np.nan
        max_relative_frag_intensity = np.nan

        compound_ref_rt_peak = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_peak
        compound_ref_rt_min = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_min
        compound_ref_rt_max = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_max
        compound_ref_mz = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz
        try:
            inchi_key = metatlas_dataset[0][compound_idx]['identification'].compound[0].inchi_key
        except:
            inchi_key = ""

        if len(msms_hits_df) == 0:
            comp_msms_hits = msms_hits_df
        else:
            comp_msms_hits = msms_hits_df[(msms_hits_df['inchi_key'] == inchi_key) \
                                          & (msms_hits_df['msms_scan'] >= compound_ref_rt_min) & (msms_hits_df['msms_scan'] <= compound_ref_rt_max) \
                                          & ((abs(msms_hits_df['measured_precursor_mz'].values.astype(float) - metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz)/metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz) \
                                             <= metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz_tolerance*1e-6)]

        for file_idx in range(len(file_names)):
            try:
                assert(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['peak_height'] > 0)
                intensities.append(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['peak_height'])
            except:# AssertionError:
                pass

            try:
                assert(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['num_ms1_datapoints'] > 0)
                rt_shifts.append(abs(compound_ref_rt_peak - metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['rt_peak']))
                mz_ppms.append(1e6*(abs(compound_ref_mz - metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['mz_centroid']) / compound_ref_mz))
            except:# AssertionError:
                pass

        if len(comp_msms_hits['score']) > 0:
            row = comp_msms_hits.loc[comp_msms_hits['score'].idxmax()]
            if np.isnan(row['msv_ref_aligned']).all():
                max_msms_score = np.nan
            else:
                max_msms_score = row['score']
            num_frag_matches = row['num_matches']

            if num_frag_matches > 1:
                msv_sample_matches = sp.partition_aligned_ms_vectors(row['msv_query_aligned'],
                                                                     row['msv_ref_aligned'])[0]
                msv_sample_matches = msv_sample_matches[:, msv_sample_matches[1].argsort()[::-1]]
                msv_sample_matches_by_intensity = msv_sample_matches[:, msv_sample_matches[1].argsort()]

                max_relative_frag_intensity = msv_sample_matches_by_intensity[1,-2] / msv_sample_matches_by_intensity[1,-1]

        try:
            max_intensity = np.nanmax(intensities)
        except ValueError:
            max_intensity = np.nan
        try:
            median_rt_shift = np.nanmedian(rt_shifts)
        except ValueError:
            median_rt_shift = np.nan
        try:
            median_mz_ppm = np.nanmedian(mz_ppms)
        except ValueError:
            median_mz_ppm = np.nan

        # assign scores
        #scores.append([metatlas_dataset[0][compound_idx]['identification'].compound[0].name,
        scores.append([compound_names[compound_idx],
                       inchi_key,
                       max_intensity,
                       median_rt_shift,
                       median_mz_ppm,
                       max_msms_score,
                       num_frag_matches,
                       max_relative_frag_intensity])

    scores_df = pd.DataFrame(scores,
                             columns=['name',
                                      'inchi_key',
                                      'max_intensity',
                                      'median_rt_shift',
                                      'median_mz_ppm',
                                      'max_msms_score',
                                      'num_frag_matches',
                                      'max_relative_frag_intensity'])

    return scores_df


def test_scores_df(scores_df,
                   min_intensity, rt_tolerance, mz_tolerance,
                   min_msms_score, allow_no_msms, min_num_frag_matches, min_relative_frag_intensity):
    """
    Returns pandas series containing boolean values for each
    compound in scores_df describing if it passes minimum requirements set by:
    'min_intensity', 'rt_tolerance','mz_tolerance', 'min_msms_score',
    'min_num_frag_matches', and 'min_relative_frag_intensity'.

    'min_intensity' <= highest intensity across all files for given compound
    'rt_tolerance' >= shift of median RT across all files for given compound to reference
    'mz_tolerance' >= ppm of median mz across all files for given compound relative to reference
    'min_msms_score' <= highest compound dot-product score across all files for given compound relative to reference
    'min_num_frag_matches' <= number of matching mzs when calculating max_msms_score
    'min_relative_frag_intensity' <= ratio of second highest to first highest intensity of matching sample mzs

    :param scores_df:

    :return test_series:
    """

    def test_row(row):
        if (row.max_intensity >= min_intensity):
            if (row.median_rt_shift <= rt_tolerance):
                if (row.median_mz_ppm <= mz_tolerance):
                    if allow_no_msms:
                        if np.isnan(row.max_msms_score):
                            return True
                    if (row.max_msms_score >= min_msms_score):
                        if (row.num_frag_matches >= min_num_frag_matches):
                            if (row.num_frag_matches <= 1):
                                return True
                            if (row.max_relative_frag_intensity >= min_relative_frag_intensity):
                                return True
        return False

    return scores_df.apply(test_row, axis=1)



def filter_atlas_and_dataset(scores_df, atlas_df, metatlas_dataset,
                             column='passing'):
    """
    Splits atlas and metatlas_dataset by compound according to if it
    passes/fails minimum requirements set by:
    'min_intensity', 'rt_tolerance','mz_tolerance', 'min_msms_score',
    'min_num_frag_matches', and 'min_relative_frag_intensity'.

    'min_intensity' <= highest intensity across all files for given compound
    'rt_tolerance' >= shift of median RT across all files for given compound to reference
    'mz_tolerance' >= ppm of median mz across all files for given compound relative to reference
    'min_msms_score' <= highest compound dot-product score across all files for given compound relative to reference
    'min_num_frag_matches' <= number of matching mzs when calculating max_msms_score
    'min_relative_frag_intensity' <= ratio of second highest to first highest intensity of matching sample mzs

    :param scores_df:
    :param atlas:
    :param metatlas_dataset:

    :return pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset:
    """

    try:
        assert(column in scores_df)
    except AssertionError:
        print 'Error: ' + column + ' not in scores_df. Either set column where pass/fail boolean values are or run test_scores_df().'
        raise
    try:
        assert(scores_df.inchi_key.tolist()
               == atlas_df.inchi_key.tolist()
               == [metatlas_dataset[0][i]['identification'].compound[0].inchi_key if len(metatlas_dataset[0][i]['identification'].compound) > 0 else ''
                                   for i in range(len(metatlas_dataset[0]))])
    except AssertionError:
        print 'Error: scores_df, atlas_df, and metatlas_dataset must have the same compounds in the same order.'
        raise

    pass_atlas_df = atlas_df[scores_df[column]]
    fail_atlas_df = atlas_df[~scores_df[column]]

    pass_dataset = []
    fail_dataset = []

    for file_idx in range(len(metatlas_dataset)):
        pass_temp = []
        fail_temp = []

        for compound_idx in range(len(metatlas_dataset[0])):
            if scores_df.iloc[compound_idx][column]:
                pass_temp.append(metatlas_dataset[file_idx][compound_idx])
            else:
                fail_temp.append(metatlas_dataset[file_idx][compound_idx])

        pass_dataset.append(pass_temp)
        fail_dataset.append(fail_temp)

    return pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset


def filter_and_output(atlas_df, metatlas_dataset, output_dir,
                      min_intensity,
                      rt_tolerance,
                      mz_tolerance,
                      min_msms_score, allow_no_msms,
                      min_num_frag_matches,  min_relative_frag_intensity,
                      num_threads=4,
                      output_pass=True, output_fail=False,
                      compress=False):

    """
    Splits atlas and metatlas_dataset by compound according to if it
    passes/fails minimum requirements set by:
    'min_intensity', 'rt_tolerance','mz_tolerance', 'min_msms_score',
    'min_num_frag_matches', and 'min_relative_frag_intensity' and
    creates error bars, chromatograms, and identification figures in output_dir.

    'min_intensity' <= highest intensity across all files for given compound
    'rt_tolerance' >= shift of median RT across all files for given compound to reference
    'mz_tolerance' >= ppm of median mz across all files for given compound relative to reference
    'min_msms_score' <= highest compound dot-product score across all files for given compound relative to reference
    'min_num_frag_matches' <= number of matching mzs when calculating max_msms_score
    'min_relative_frag_intensity' <= ratio of second highest to first highest intensity of matching sample mzs
    'num_threads' = number of threads to use in multiprocessing

    Returns the unfiltered metatlas dataset and filtered dataset that can be used for downstream processing steps.

    :param atlas:
    :param groups:
    :param output_dir:
    """

    with open(os.path.join(output_dir, 'test_parameters.txt'), 'w') as f:
        f.write('min_intensity=' + str(min_intensity) + '\n' +
                'rt_tolerance=' + str(rt_tolerance) + '\n' +
                'mz_tolerance=' + str(mz_tolerance) + '\n' +
                'min_msms_score=' + str(min_msms_score) + '\n' +
                'allow_no_msms=' + str(allow_no_msms) + '\n' +
                'min_num_frag_matches=' + str(min_num_frag_matches) + '\n' +
                'min_relative_frag_intensity=' + str(min_relative_frag_intensity))

    print 'making scores_df'
    # scores compounds in metatlas dataset
    scores_df = make_scores_df(metatlas_dataset)

    print 'testing and making compound_scores.csv'
    # scores dataframe
    scores_df['passing'] = test_scores_df(scores_df,
                                          min_intensity, rt_tolerance, mz_tolerance,
                                          min_msms_score, allow_no_msms, min_num_frag_matches, min_relative_frag_intensity)
    scores_df.to_csv(os.path.join(output_dir, 'stats_tables', 'compound_scores.csv'))

    print 'filtering atlas and dataset'
    # filter dataset by scores
    pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset = filter_atlas_and_dataset(scores_df, atlas_df, metatlas_dataset)

    outputs = []

    if output_pass:
        try:
            pass_dataset[0][0]['data']
            outputs.append((pass_atlas_df, pass_dataset, os.path.join(output_dir, 'pass')))
        except:
            pass

    if output_fail:
        try:
            fail_dataset[0][0]['data']
            outputs.append((fail_atlas_df, fail_dataset, os.path.join(output_dir, 'fail')))
        except:
            pass

    for atlas_df, filtered_dataset, output_dir in outputs:

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        print 'saving atlas'
        atlas_df.to_csv(os.path.join(output_dir, 'filtered_atlas_export.csv'))

        print 'making info tables'
        peak_height = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_height' , output_loc=os.path.join(output_dir,'sheets'))
        peak_area = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_area' , output_loc=os.path.join(output_dir,'sheets'))
        mz_peak = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_peak' , output_loc=os.path.join(output_dir,'sheets'))
        rt_peak = dp.make_output_dataframe(input_fname = '', input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [],fieldname='rt_peak' , output_loc=os.path.join(output_dir,'sheets'))
        mz_centroid = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_centroid' , output_loc=os.path.join(output_dir,'sheets'))
        rt_centroid = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='rt_peak' , output_loc=os.path.join(output_dir,'sheets'))

        print 'making error bars'
        #Error bars
        peak_height = dp.make_output_dataframe(input_fname='', input_dataset=filtered_dataset, include_lcmsruns=[], exclude_lcmsruns=[], fieldname='peak_height')
        dp.plot_errorbar_plots(peak_height, output_loc=os.path.join(output_dir, 'error_bar_peak_height'))

        print 'making identification figures'
        #Identification figures
        dp.make_identification_figure_v2(input_dataset=filtered_dataset, include_lcmsruns = [],exclude_lcmsruns=[], output_loc=os.path.join(output_dir, 'identification'))

        print 'making chromatograms'
        # Chromatograms
        group = 'sort'  # 'page' or 'index' or 'sort' or None
        save = True
        share_y = True

        file_names = ma_data.get_file_names(filtered_dataset)
        compound_names = ma_data.get_compound_names(filtered_dataset)[0]
        args_list = []

        chromatogram_str = 'compound_chromatograms'

        if not os.path.exists(os.path.join(output_dir, chromatogram_str)):
            os.makedirs(os.path.join(output_dir, chromatogram_str))

        for compound_idx, my_compound in enumerate(compound_names):
            my_data = list()
            for file_idx, my_file in enumerate(file_names):
                my_data.append(filtered_dataset[file_idx][compound_idx])
            kwargs = {'data': my_data,
                      'file_name': os.path.join(output_dir, chromatogram_str, my_compound+'.pdf'),
                      'group': group,
                      'save': save,
                      'share_y': share_y,
                      'names': file_names}
            args_list.append(kwargs)

        pool = mp.Pool(processes=min(num_threads, len(filtered_dataset[0])))
        pool.map(cpp.chromplotplus, args_list)
        pool.close()
        pool.terminate()

    print 'done'
    return pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset
