#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


from .issue_generator_clang_tool import ClangToolIssueGenerator
import copy
import difflib


class ClangFormatIssueGenerator(ClangToolIssueGenerator):
    """issue generator for clang format"""

    def __init__(self, **kwargs):
        kwargs["clang_tool_executable"] = kwargs.pop(
            "clang_format_executable", "clang-format"
        )
        super().__init__("clang-format", **kwargs)

    def _get_classifiaction(self, issue):
        return (["Style"], "info")

    def _generate_issues_for_file(self, filename):
        with open(filename, "r") as fh:
            unformated_file = fh.readlines()

        command = copy.copy(self._command)
        command.append(filename)
        if self._current_file is not None:
            self._current_file.set_description_str(f"Current file: {filename}")

        formated_file = list(self._run_command(command))

        issues = []
        s = difflib.SequenceMatcher(None, unformated_file, formated_file)
        for tag, i1, i2, j1, j2 in s.get_opcodes():
            description = ""
            if tag == "equal":
                continue
            if tag == "delete":
                description = f"```suggestion:-0+{i2 - i1 - 1}\n```"
            if tag == "insert":
                description = f"""```suggestion:-0+0\n{"".join(unformated_file[i1:i2])}{"".join(formated_file[j1:j2])}```"""  # noqa
            if tag == "replace":
                description = f"""```suggestion:-0+{i2 - i1 - 1}\n{"".join(formated_file[j1:j2])}```"""  # noqa

            issue = {
                "name": f"""clang-format:{tag}""",
                "description": "".join(description),
                "file": filename,
                "line": i1 + 1,  # lines start at 1 not 0
                "column": 1,
                "end_line": i2 + 1,
            }
            issues.append(issue)
        return issues
