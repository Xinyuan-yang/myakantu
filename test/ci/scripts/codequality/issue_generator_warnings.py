#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


from . import print_info
from .issue_generator import IssueGenerator
import re
import warning_parser as warn


class WarningsIssueGenerator(IssueGenerator):
    """Main class to run and convert the results of clang-tidy to the
    code-quality format.
    """

    CLASSIFICATIONS = {
        "gcc": {
            "uninitialized": {
                "categories": ["Bug Risk"],
                "severity": "major",
            },
            "sign-compare": {"categories": ["Bug Risk"], "severity": "minor"},
        },
    }

    def __init__(self, **kwargs):
        """Construct the class."""
        super().__init__(**kwargs)
        files = kwargs.pop("files")
        compiler = kwargs.pop("compiler", "auto")

        self._input_files = {}
        if compiler == "auto":
            compiler_re = re.compile(r".*build.*(gcc|clang).*-(out|err)\.log")

            for _file in files:
                match = compiler_re.search(_file)
                if match:
                    compiler = match.group(1)
                    if compiler not in self._input_files:
                        self._input_files[compiler] = []

                    self._input_files[compiler].append(_file)
                else:
                    print_info(f"Skipped {_file}, could not determine compiler")
        else:
            self._input_files[compiler] = files

    def generate_issues(self):
        """Parse warning files."""
        for compiler, _files in self._input_files.items():
            for _file in _files:
                warnings = warn.get_warnings(_file, compiler)
                for warning in warnings:
                    issue = {
                        "name": f"warning:{compiler}-{warning.get_category()}",
                        "description": warning.get_message(),
                        "file": warning.get_filepath(),
                        "line": warning.get_line(),
                        "column": warning.get_column(),
                        "raw": warning,
                    }

                    self.add_issue(issue)

    def _get_classifiaction(self, warning):
        categories = ["Clarity"]
        severity = "info"

        tool = warning["raw"].get_tool()
        cat = warning["raw"].get_category()

        if tool in self.CLASSIFICATIONS:
            classifications = self.CLASSIFICATIONS[tool]
            if cat in classifications:
                categories = classifications[cat]["categories"]
                severity = classifications[cat]["severity"]

        return (categories, severity)
