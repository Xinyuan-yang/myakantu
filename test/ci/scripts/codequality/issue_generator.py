#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


from . import print_debug, print_info
import hashlib
import os
import re
import copy


class IssueGenerator:
    """Interface for the issue generators"""

    def __init__(self, **kwargs):
        self._files = kwargs.pop("file_list", [])

        excludes = kwargs.pop("excludes", None)
        if excludes is None:
            excludes = []

        extensions = kwargs.pop("extensions", None)
        if extensions is None:
            extensions = [".cc", ".hh"]

        self._extensions = [
            re.compile(r"\{}$".format(extension)) for extension in extensions
        ]

        self._exclude_patterns = [re.compile(exclude) for exclude in excludes]

        self._issues = {}
        self._filter_file_list()

    def _filter_file_list(self):
        file_list = copy.copy(self._files)
        self._files = []
        for filename in file_list:
            filename = os.path.relpath(filename)
            need_exclude = self._need_exclude(filename)
            if need_exclude:
                print_debug(f"exluding file: {filename}")
                continue
            print_info(f"adding file: {filename}")
            self._files.append(filename)

    def _need_exclude(self, filename):
        need_exclude = False
        for pattern in self._exclude_patterns:
            match = pattern.search(filename)
            need_exclude |= bool(match)

        match_extension = False
        for extension in self._extensions:
            match = extension.search(filename)
            match_extension |= bool(match)

        need_exclude |= not match_extension
        return need_exclude

    def add_issue(self, unfmt_issue):
        """add an issue to the list if not already present"""
        issue = self._format_issue(unfmt_issue)

        filepath = issue["location"]["path"]
        if self._need_exclude(filepath):
            return

        if issue["fingerprint"] in self._issues:
            return

        self._issues[issue["fingerprint"]] = issue

    @property
    def issues(self):
        """get the list of registered issues"""
        return list(self._issues.values())

    def _format_issue(self, unfmt_issue):
        filepath = os.path.relpath(unfmt_issue["file"])
        issue = {
            "type": "issue",
            "check_name": unfmt_issue["name"],
            "description": (
                f"""[{unfmt_issue["name"]}] {unfmt_issue["description"]}"""
            ),
            "location": {
                "path": filepath,
                "lines": {
                    "begin": unfmt_issue["line"],
                    "end": unfmt_issue["line"],
                },
                "positions": {
                    "begin": {
                        "line": unfmt_issue["line"],
                        "column": unfmt_issue["column"],
                    },
                    "end": {
                        "line": unfmt_issue["line"],
                        "column": unfmt_issue["column"],
                    },
                },
            },
        }

        if "end_line" in unfmt_issue:
            issue["location"]["positions"]["end"] = {
                "line": unfmt_issue["end_line"],
                "column": unfmt_issue["column"],
            }
            issue["location"]["lines"]["end"] = unfmt_issue["end_line"]

        issue["fingerprint"] = hashlib.md5(
            "{file}:{line}:{column}:{type}".format(
                file=filepath,
                line=unfmt_issue["line"],
                column=unfmt_issue["column"],
                type=unfmt_issue["name"],
            ).encode()
        ).hexdigest()

        issue["categories"], issue["severity"] = self._get_classifiaction(
            unfmt_issue
        )

        print_debug(issue)
        return issue
