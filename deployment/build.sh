#!/bin/bash
# Install distutils
python3 -m ensurepip --upgrade
python3 -m pip install --upgrade pip setuptools wheel

# Install dependencies
python3 -m pip install -r requirements.txt