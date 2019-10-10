import smtplib
import sys

from datetime import datetime, time as dtime

def send_mail(subject, username, body, force=False):
    """Send the mail only once per day."""

    now = datetime.now()
    if force or dtime(00, 00) <= now.time() <= dtime(00, 10):
        sender = 'pasteur@nersc.gov'
        receivers = ['%s@nersc.gov' % username, '%s@nersc.gov' % 'bpb']
        message = """\
From: %s
To: %s
Subject: %s

%s
        """ % (sender, ", ".join(receivers), subject, body)
        try:
            smtpObj = smtplib.SMTP('localhost')
            smtpObj.sendmail(sender, receivers, message)
            sys.stdout.write("Successfully sent email to %s\n" % username)
            sys.stdout.flush()
        except smtplib.SMTPException:
            sys.stderr.write("Error: unable to send email to %s\n" % username)
            sys.stdout.flush()
