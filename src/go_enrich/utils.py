#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
General utility functions
"""

import re
import os
import sys


def parse_obo_file(obo_path):
    """
    Parse the go-basic.obo file and return a dict {GO_id: name}.
    Only extracts 'id' and 'name' fields, ignoring obsolete terms.
    """
    go_map = {}
    with open(obo_path, 'r', encoding='utf-8') as f:
        current_id = None
        current_name = None
        is_obsolete = False
        for line in f:
            line = line.strip()
            if not line:
                # Empty line marks the end of a Term
                if current_id and current_name and not is_obsolete:
                    go_map[current_id] = current_name
                current_id = None
                current_name = None
                is_obsolete = False
                continue
            if line.startswith('id: GO:'):
                current_id = line.split('id: ')[1]
            elif line.startswith('name: '):
                current_name = line.split('name: ')[1]
            elif line.startswith('is_obsolete: true'):
                is_obsolete = True
        # Handle case where file ends without a blank line
        if current_id and current_name and not is_obsolete:
            go_map[current_id] = current_name
    return go_map


def get_default_obo_path():
    """
    Return the path to the built-in go-basic.obo file inside the package data directory.
    Handles both normal development and PyInstaller bundled environment.
    Returns None if the file does not exist.
    """
    # Check if running as a PyInstaller bundle
    if getattr(sys, 'frozen', False):
        # Running in a PyInstaller bundle
        base_path = sys._MEIPASS
    else:
        # Running in normal Python environment
        base_path = os.path.dirname(os.path.abspath(__file__))

    obo_path = os.path.join(base_path, 'data', 'go-basic.obo')
    if os.path.exists(obo_path):
        return obo_path
    return None