# app.py
import os

from flask import Flask

APP_ROOT = os.path.dirname(os.path.realpath(__file__))
DEBUG = False

class CustomFlask(Flask):
  jinja_options = Flask.jinja_options.copy()
  jinja_options.update(dict(
    block_start_string='(%',
    block_end_string='%)',
    variable_start_string='((',
    variable_end_string='))',
    comment_start_string='(#',
    comment_end_string='#)',
  ))

app = CustomFlask(__name__)
app.config.from_object(__name__)
UPLOAD_FOLDER = './tempuploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

try:
    os.mkdir(UPLOAD_FOLDER)
except:
    print("cannot make folder", UPLOAD_FOLDER)

try:
    os.mkdir(os.path.join(UPLOAD_FOLDER, "reference_spectra"))
except:
    print("cannot make folder reference_spectra")
