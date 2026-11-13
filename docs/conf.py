#
# Configuration file for the Sphinx documentation builder.

from datetime import datetime
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = "NBSOS"
author = "The World Bank, Global Program on Nature-Based Solutions for Climate Resilience"
version = "0.0.1"
release = version

now = datetime.now()
today = f"{now.year}-{now.month:02}-{now.day:02} {now.hour:02}H{now.minute:02}"
copyright = f"{now.year}, {author}"

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

numfig = True

numfig_format = {
    'figure': 'Figure %s.'
    }

numfig_secnum_depth = 1

pdf_documents = [('index', u'Urban_NBSOS', u'URBAN NBSOS', u'GPNBS GFDRR WB'),]

  # index - master document
  # rst2pdf - name of the generated pdf
  # Sample rst2pdf doc - title of the pdf
  # Your Name - author name in the pdf

master_doc = "index"
language = "en"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".venv"]
pygments_style = "sphinx"
extensions = [
    "sphinx.ext.intersphinx", 
    "sphinx.ext.autodoc",
    "sphinx_rtd_theme",
    "rst2pdf.pdfbuilder",
    "docxbuilder"
]

todo_include_todos = False
html_theme = "sphinx_rtd_theme"
