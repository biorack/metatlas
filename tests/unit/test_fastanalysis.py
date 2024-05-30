'''unit tests of fastanalysis functions'''

from metatlas.tools import fastanalysis
import numpy as np
import pandas as pd
import statistics

def test_calculate_compound_total_score001():
    
    # Create a mock DataFrame
    final_df = pd.DataFrame({
        'total_score': [0],
        'msms_quality': [0],
        'msi_level': ['']
    })

    # Define the compound index and quality scores
    compound_idx = 0
    quality_scores = [1, 1, 1]

    # Call the function with the test parameters
    result_df = fastanalysis.calculate_compound_total_score(final_df, compound_idx, quality_scores)

    # Assert that the total_score is the sum of the quality scores
    assert result_df.loc[compound_idx, 'msi_level'] == "Exceeds Level 1"


def test_calculate_compound_total_score002():
    
    # Create a mock DataFrame
    final_df = pd.DataFrame({
        'total_score': [0],
        'msms_quality': [0],
        'msi_level': ['']
    })

    # Define the compound index and quality scores
    compound_idx = 0
    quality_scores = [0, 1, 1]

    # Call the function with the test parameters
    result_df = fastanalysis.calculate_compound_total_score(final_df, compound_idx, quality_scores)

    # Assert that the total_score is the sum of the quality scores
    assert result_df.loc[compound_idx, 'msi_level'] == "Level 1"


def test_calculate_compound_total_score003():
    
    # Create a mock DataFrame
    final_df = pd.DataFrame({
        'total_score': [0],
        'msms_quality': [0],
        'msi_level': ['']
    })

    # Define the compound index and quality scores
    compound_idx = 0
    quality_scores = [0, 0, 0]

    # Call the function with the test parameters
    result_df = fastanalysis.calculate_compound_total_score(final_df, compound_idx, quality_scores)

    # Assert that the total_score is the sum of the quality scores
    assert result_df.loc[compound_idx, 'msi_level'] == "putative"


def test_calculate_compound_total_score004():
    
    # Create a mock DataFrame
    final_df = pd.DataFrame({
        'total_score': [0],
        'msms_quality': [-1],
        'msi_level': ['']
    })

    # Define the compound index and quality scores
    compound_idx = 0
    quality_scores = [-1, 0, 0]

    # Call the function with the test parameters
    result_df = fastanalysis.calculate_compound_total_score(final_df, compound_idx, quality_scores)

    # Assert that the total_score is the sum of the quality scores
    assert result_df.loc[compound_idx, 'msi_level'] == "REMOVE, INVALIDATED BY BAD MSMS MATCH"