""" stand alone utility functions """


def or_default(none_or_value, default):
    """
    inputs:
        none_or_value: variable to test
        default: value to return if none_or_value is None
    """
    return none_or_value if none_or_value is not None else default
