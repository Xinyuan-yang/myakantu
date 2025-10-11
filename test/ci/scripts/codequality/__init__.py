#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import json as _cq_json
import sys as _cq_sys

try:
    from termcolor import colored as _cq_colored
except ImportError:

    def _cq_colored(text, color):  # pylint: disable=unused-argument
        """fallback function for termcolor.colored"""
        return text


def print_debug(message):
    """helper function to print debug messages"""
    print(
        f"Debug: {_cq_colored(message, 'red')}",
        file=_cq_sys.stderr,
        flush=True,
    )


def print_info(message):
    """helper function to print info messages"""
    print(
        f"Info: {_cq_colored(message, 'blue')}",
        file=_cq_sys.stderr,
        flush=True,
    )


def run(cmd, **kwargs):
    from .issue_generator_clang_tidy import ClangTidyIssueGenerator  # NOQA
    from .issue_generator_clang_format import ClangFormatIssueGenerator  # NOQA
    from .issue_generator_warnings import WarningsIssueGenerator  # NOQA

    if cmd == "clang_tidy":
        tool = ClangTidyIssueGenerator(**kwargs)
    elif cmd == "clang_format":
        tool = ClangFormatIssueGenerator(**kwargs)
    elif cmd == "warnings":
        tool = WarningsIssueGenerator(**kwargs)

    tool.generate_issues()

    print(_cq_json.dumps(tool.issues))
