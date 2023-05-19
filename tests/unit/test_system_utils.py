""" unit testing of system_utils functions """
# pylint: disable=missing-function-docstring

import pytest

from metatlas.io import system_utils


def test_send_mail01():
    with pytest.raises(ValueError):
        system_utils.send_mail("subject", "-fake_switch not_an_email_address", "body")


def test_send_mail02():
    with pytest.raises(ValueError):
        system_utils.send_mail("subject", "not_an_email_address", "body")


def test_send_mail03():
    with pytest.raises(ValueError):
        system_utils.send_mail("subject", [], "body")


def test_send_mail04():
    with pytest.raises(ValueError):
        system_utils.send_mail("subject", "", "body")
