""" Object Oriented Groups"""
import logging
import pandas as pd

from metatlas.datastructures import metatlas_objects as metob
from metatlas.plots import dill2plots as dp

logger = logging.getLogger(__name__)


class Groups:
    """Groups of LCMS run to define related experimental samples"""
    def __init__(self, experiment, analysis_id, controlled_vocab, exclude_files):
        self.experiment = experiment
        self.analysis_id = analysis_id
        self.controlled_vocab = controlled_vocab
        self.exclude_files = exclude_files

    def query_lcmsruns(self, most_recent=True):
        """Get LCMS runs from DB matching experiment"""
        return dp.get_metatlas_files(experiment=self.experiment, name='%', most_recent=most_recent)

    def get_lcmsruns_df(self, most_recent=True):
        """Returns a pandas DataFrame with lcmsrun matching self.experiment"""
        files = dp.get_metatlas_files(experiment=self.experiment, name='%', most_recent=most_recent)
        logger.info("Number of LCMS output files matching '%s' is: %d.", self.experiment, len(files))
        return metob.to_dataframe(files)

    lcmsruns_df = property(get_lcmsruns_df)

    def get_lcmsruns_short_names(self, fields=None):
        """
        Querys DB for lcms filenames from self.experiment and returns
        a pandas DataFrame containing identifiers for each file
        inputs:
            fields: optional dict with column names as key
                    and list of lcms filename metadata fields positions as value
        """
        if fields is None:
            fields = {'full_filename': range(16),
                      'sample_treatment': [12],
                      'short_filename': [0, 2, 4, 5, 7, 9, 14],
                      'short_samplename': [0, 2, 4, 5, 7, 9, 14],
                      }
        out = pd.DataFrame(columns=fields.keys())
        for i, lcms_file in enumerate(self.query_lcmsruns()):
            tokens = lcms_file.name.split('.')[0].split('_')
            for name, idxs in fields.items():
                out.loc[i, name] = "_".join([tokens[n] for n in idxs])
            out.loc[i, 'last_modified'] = pd.to_datetime(lcms_file.last_modified, unit='s')
        out.sort_values(by='last_modified', inplace=True)
        out.drop(columns=['last_modified'], inplace=True)
        out.drop_duplicates(subset=['full_filename'], keep='last', inplace=True)
        out.set_index('full_filename', inplace=True)
        return out

    lcmsruns_short_names = property(get_lcmsruns_short_names)

    def get_group_name(self, base_filename):
        """Returns dict with keys group and short_name corresponding to base_filename"""
        indices = [i for i, s in enumerate(self.controlled_vocab) if s.lower() in base_filename.lower()]
        tokens = base_filename.split('_')
        prefix = '_'.join(tokens[:11])
        suffix = self.controlled_vocab[indices[0]].lstrip('_') if indices else tokens[12]
        group_name = f"{prefix}_{self.analysis_id}_{suffix}"
        short_name = f"{tokens[9]}_{suffix}"  # Prepending POL to short_name
        return {'group': group_name, 'short_name': short_name}

    @property
    def _get_files_dict(self):
        """
        Queries DB for all lcmsruns matching the class properties.
        Returns a dict of dicts where keys are filenames minus extensions and values are
        dicts with keys: name (filename with extension), group, and short_name
        """
        file_dict = {}
        for lcms_file in self.query_lcmsruns():
            if not any(map(lcms_file.name.__contains__, self.exclude_files)):
                base_name = lcms_file.name.split('.')[0]
                file_dict[base_name] = {'name': lcms_file.name, **self.get_group_name(base_name)}
        return file_dict

    @property
    def df(self):  # pylint: disable=invalid-name
        """Returns pandas Dataframe with one group per row"""
        out = pd.DataFrame(self._get_files_dict).T
        out.index.name = 'filename'
        out.reset_index(inplace=True)
        out.drop(columns=['name'], inplace=True)
        return out

    @property
    def group_objects(self):
        """Returns a list of Group objects"""
        file_dict = self._get_files_dict
        out = []
        for group_name, values in self.df.to_dict('index'):
            out.append(metob.Group(name=group_name,
                                   short_name=values['short_name'],
                                   items=[file_value['name']
                                          for file_value in file_dict.values()
                                          if file_value['group'] == group_name]))
        return out
