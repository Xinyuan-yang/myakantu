#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2022-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import os
import re
import sys
import subprocess


def run_git_command(args):
    """Run git commands and capture outputs."""
    git_dir = os.path.realpath(
        os.path.join(os.path.dirname(__file__), os.pardir)
    )

    cmd = ["git"] + args
    p = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=git_dir
    )
    stdout = p.communicate()[0].strip().decode()
    if p.returncode != 0:
        return None, p.returncode
    return stdout, p.returncode


def _split_git_describe(describe):
    describe_mo = re.search(
        r"^(?P<tag>.+)"
        r"-(?P<distance>\d+)"
        r"-g(?P<short>[0-9a-f]+)"
        r"(-(?P<dirty>dirty))?$",
        describe,
    )

    if describe_mo:
        pieces = {}
        pieces["tag"] = describe_mo.group("tag")
        # distance: number of commits since tag
        pieces["distance"] = int(describe_mo.group("distance"))

        # commit: short hex revision ID
        pieces["short"] = describe_mo.group("short")
        if describe_mo.group("dirty"):
            pieces["dirty"] = describe_mo.group("dirty")

        return pieces

    return None


def _eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


semver_re = re.compile(
    r"^(?P<major>0|[1-9]\d*)"
    r"(\.(?P<minor>0|[1-9]\d*))?"
    r"(\.(?P<patch>0|[1-9]\d*))?"
    r"(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)"
    r"(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?"
    r"(?:\+(?P<build>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"
)


def _parse_semver(version):
    semver_mo = semver_re.search(version)
    if not semver_mo:
        return {}

    pieces = {}
    for p in ["major", "minor", "patch", "prerelease", "build"]:
        if semver_mo.group(p):
            pieces[p] = semver_mo.group(p)
    return pieces


def get_git_version():
    """Get the version from the git repository."""
    out, rc = run_git_command(["rev-parse", "--git-dir"])
    if rc != 0:
        return None

    git_describe, rc = run_git_command(
        ["describe", "--tags", "--dirty", "--always", "--match", "v*"]
    )

    _eprint(f"git describe {git_describe}")
    if "-g" in git_describe:
        # TAG-DISTANCE-gHEX
        pieces = _split_git_describe(git_describe)
    else:
        # tag only or no tag and hash
        pieces = {"tag": git_describe}

    # major.minor.patch-prerelease+build
    if not pieces or ("tag" not in pieces):
        return None

    semver_mo = semver_re.search(pieces["tag"][1:])
    if not semver_mo:
        return pieces

    for p in ["major", "minor", "patch", "prerelease", "build"]:
        if semver_mo.group(p):
            pieces[p] = semver_mo.group(p)

    return pieces


def get_git_attributes_version():
    """Get the version from the attributes set a `git archive`."""
    file_dir = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
    attributes = None
    pieces = None

    with open(os.path.join(file_dir, "git_info")) as fh:
        describe_re = re.compile(r"describe: ([^$].*[^$])")
        for line in fh:
            mo = describe_re.search(line)
            if mo:
                attributes = mo.group(1)
                break

    if attributes:
        pieces = _split_git_describe(attributes)

    return pieces


def get_ci_version():
    """Get extra information from CI context."""
    pieces = None
    if "CI_AKANTU_INSTALL_PREFIX" not in os.environ:
        return None

    ci_akantu_install_prefix = os.environ["CI_AKANTU_INSTALL_PREFIX"]
    akantu_dir = os.path.join(
        ci_akantu_install_prefix, "lib", "cmake", "Akantu"
    )
    cmake_config = os.path.join(akantu_dir, "AkantuConfig.cmake")

    if not os.path.exists(cmake_config):
        return None

    version = None
    with open(cmake_config, "r") as fh:
        version_re = re.compile(r"^set\(AKANTU_VERSION (.*)\)$")
        for line in fh:
            version_mo = version_re.search(line)
            if version_mo:
                version = version_mo.group(1)
                break

    if not version:
        return None

    pieces = _parse_semver(version)
    return pieces


def get_mr_info():
    pieces = {}
    if "CI_MERGE_REQUEST_IID" in os.environ:
        pieces = {"build": "mr" + os.environ["CI_MERGE_REQUEST_IID"]}
    return pieces


def get_version_file():
    """Get the version directly from the VERSION file."""
    version_path = os.path.join(
        os.path.realpath(
            os.path.abspath(
                os.path.join(os.path.dirname(__file__), os.path.pardir)
            )
        ),
        "VERSION",
    )

    if not os.path.exists(version_path):
        return None

    version = None
    with open(version_path, "r") as fh:
        version = fh.readline()

    if not version:
        return None

    pieces = _parse_semver(version)
    return pieces


def get_version():
    """Combine all the version determination functions."""
    pieces = None

    if not pieces:
        pieces = get_version_file()
        _eprint(f"pieces from version file: {pieces}")

    if not pieces:
        pieces = get_ci_version()
        _eprint(f"pieces from ci_version: {pieces}")

    if not pieces:
        pieces = get_git_version()
        _eprint(f"pieces from git_version: {pieces}")

    if not pieces:
        pieces = get_git_attributes_version()
        _eprint(f"pieces from git attributes: {pieces}")

    if not pieces:
        raise Exception("No version could be determined")

    semver_build = []
    if "build" in pieces:
        semver_build.append(pieces["build"])

    ci_build_pieces = get_mr_info()
    if "build" in ci_build_pieces:
        _eprint(f"pieces from version MR: {ci_build_pieces}")
        semver_build.append(ci_build_pieces["build"])

    if "distance" in pieces:
        semver_build.extend([str(pieces["distance"]), "g" + pieces["short"]])
        if "dirty" in pieces and pieces["dirty"]:
            semver_build.append(pieces["dirty"])

    if semver_build:
        pieces["build_part"] = "+" + ".".join(semver_build)
    else:
        pieces["build_part"] = ""

    if "prerelease" in pieces:
        pieces["prerelease"] = "-" + pieces["prerelease"]
    else:
        pieces["prerelease"] = ""

    semver = "{major}.{minor}.{patch}{prerelease}{build_part}".format(**pieces)

    return semver


if __name__ == "__main__":
    print(get_version())
