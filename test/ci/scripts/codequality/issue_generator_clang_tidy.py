#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


from . import print_debug, print_info
from .issue_generator_clang_tool import ClangToolIssueGenerator
import os
import re
import copy
import json
import subprocess
import multiprocessing as mp


class ClangTidyIssueGenerator(ClangToolIssueGenerator):
    """issue generator for clang tidy"""

    # 7-bit C1 ANSI sequences
    ANSI_ESCAPE = re.compile(
        r"""
        \x1B  # ESC
        (?:   # 7-bit C1 Fe (except CSI)
            [@-Z\\-_]
        |     # or [ for CSI, followed by a control sequence
            \[
            [0-?]*  # Parameter bytes
            [ -/]*  # Intermediate bytes
            [@-~]   # Final byte
        )
    """,
        re.VERBOSE,
    )

    ISSUE_PARSE = re.compile(
        r"(?P<file>.*\.(cc|hh)):(?P<line>[0-9]+):(?P<column>[0-9]+): (warning|error): (?P<description>.*) \[(?P<name>.*)\]"
    )  # NOQA pylint: disable=line-too-long

    CLASSIFICATIONS = {
        "bugprone": {
            "categories": ["Bug Risk"],
            "severity": "major",
        },
        "modernize": {
            "categories": ["Clarity", "Compatibility", "Style"],
            "severity": "info",
        },
        "mpi": {
            "categories": ["Bug Risk", "Performance"],
            "severity": "critical",
        },
        "openmp": {
            "categories": ["Bug Risk", "Performance"],
            "severity": "critical",
        },
        "performance": {
            "categories": ["Performance"],
            "severity": "minor",
        },
        "readability": {
            "categories": ["Clarity", "Style"],
            "severity": "info",
        },
    }

    def __init__(self, **kwargs):
        kwargs["clang_tool_executable"] = kwargs.pop(
            "clang_tidy_executable", "clang-tidy"
        )
        super().__init__("clang-tidy", need_compiledb=True, **kwargs)
        self._num_threads = min(2, mp.cpu_count())

    def _get_classifiaction(self, issue):
        type_ = issue["type"]
        categories = ["Bug Risk"]
        severity = "blocker"

        if type_ in self.CLASSIFICATIONS:
            categories = self.CLASSIFICATIONS[type_]["categories"]
            severity = self.CLASSIFICATIONS[type_]["severity"]
        elif type_[0] == "clang":
            if type_[1] == "diagnostic":
                categories = ["Bug Risk"]
                severity = "blocker"
            elif type_[1] == "analyzer":
                categories = ["Bug Risk"]
                severity = "major"

        return (categories, severity)

    def _generate_issues_for_file(self, filename):
        command = copy.copy(self._command)
        command.append(filename)
        if self._current_file is not None:
            self._current_file.set_description_str(f"Current file: {filename}")

        issues = []
        issue = {}
        for line in self._run_command(command):
            line = line.rstrip()
            match = self.ISSUE_PARSE.match(line)
            if match:
                if len(issue) != 0:
                    issues.append(issue)
                issue = match.groupdict()
                issue["type"] = issue["name"]
                issue["name"] = f"""clang-tidy:{issue["name"]}"""

                # print_debug(f'[clang-tidy] new issue: {line}')
            elif issue:
                if "content" in issue:
                    issue["content"].append(line)
                    # print_debug(f'[clang-tidy] more extra content: {line}')
                else:
                    issue["content"] = [line]
                    # print_debug(f'[clang-tidy] extra content: {line}')
            if len(issue) != 0:
                issues.append(issue)
        return issues
