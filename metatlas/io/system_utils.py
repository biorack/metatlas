import getpass
import re
import subprocess

from typing import Optional, Sequence, Union

EMAIL_RE = re.compile(r"[^@\s]+@[^@\s]+\.[a-zA-Z0-9]+$")


def is_bad_email_address(address: str) -> bool:
    """
    Returns True only if address is obviously bad.
    Makes sure a command line switch ('-foo') isn't passed of as an email address
    """
    clean = address.strip()
    return clean.startswith("-") or not EMAIL_RE.match(clean)


def send_mail(
    subject: str, recipients: Union[str, Sequence[str]], body: str, sender: Optional[str] = None
) -> None:
    """Sends an email using the mailx command"""
    sender = getpass.getuser() if sender is None else sender
    recipients = [recipients] if isinstance(recipients, str) else recipients
    if recipients is None or len(recipients) == 0:
        raise ValueError("No recipients supplied")
    for address in recipients:
        if is_bad_email_address(address):
            raise ValueError(f"'{address}' does not appear to be an email address")
    subprocess.run(["mailx", "-s", subject] + list(recipients), text=True, input=body, check=True)
