"""
This is the access point to the parser for the scala layer in the
nomad project.
"""
from __future__ import absolute_import
import sys
import setup_paths

from LammpsParser import LammpsParser
from nomadcore.parser_backend import JsonParseEventsWriterBackend


if __name__ == "__main__":

    # Initialise the parser with the main filename and a JSON backend
    main_file = sys.argv[1]
    parser = LammpsParser(main_file, backend=JsonParseEventsWriterBackend)
    parser.parse()
