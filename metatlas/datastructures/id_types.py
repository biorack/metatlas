""" types for use with MetatlasDataset """

from typing import List, NewType, Optional, TypedDict

from traitlets import HasTraits, TraitType

import metatlas.datastructures.metatlas_objects as metob

GroupList = Optional[List[metob.Group]]
LcmsRunsList = Optional[List[metob.LcmsRun]]
FileMatchList = List[str]
GroupMatchList = List[str]

Polarity = NewType("Polarity", str)
ShortPolarity = NewType("ShortPolarity", str)
Experiment = NewType("Experiment", str)
OutputType = NewType("OutputType", str)
AnalysisNumber = NewType("AnalysisNumber", int)
PathString = NewType("PathString", str)

OUTPUT_TYPES = [
    OutputType("ISTDsEtc"),
    OutputType("FinalEMA-HILIC"),
    OutputType("FinalEMA-C18"),
    OutputType("data_QC"),
    OutputType("other"),
]
POLARITIES = [Polarity("positive"), Polarity("negative"), Polarity("fast-polarity-switching")]
SHORT_POLARITIES = {
    Polarity("positive"): ShortPolarity("POS"),
    Polarity("negative"): ShortPolarity("NEG"),
    Polarity("fast-polarity-switching"): ShortPolarity("FPS"),
}


class Proposal(TypedDict):
    """for use with traitlets.validate"""

    owner: HasTraits
    value: object
    trait: TraitType


class LcmsRunDict(TypedDict):
    """part of return type for AnalysisIdentifiers._files_dict"""

    object: metob.LcmsRun
    group: str
    short_name: str