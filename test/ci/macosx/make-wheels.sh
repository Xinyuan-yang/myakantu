#!/usr/bin/env bash
set -eo pipefail

set +x

source spack/share/spack/setup-env.sh
spack env activate macosx

python3 -m venv venv
source venv/bin/activate

pip3 install -r dev-requirements.txt

pip3 wheel . --no-deps -w dist/
delocate-wheel -w wheelhouse_macos/ --ignore-missing-dependencies -v dist/*.whl
