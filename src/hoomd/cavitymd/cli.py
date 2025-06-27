#!/usr/bin/env python3
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Command-line interface for cavity molecular dynamics simulations."""

import sys
from .experiments import run_cavity_experiments


def main():
    """Main entry point for cavity MD experiments."""
    return run_cavity_experiments()


if __name__ == '__main__':
    sys.exit(main()) 