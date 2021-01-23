from os.path import (join, abspath, dirname)
import json
import subprocess
from glob import glob

from cookiecutter.main import cookiecutter

HERE = abspath(join(dirname(__file__)))
BUILD = join(HERE, "build")

cookiecutter_path = join(HERE, "cookiecutter")

with open(join(cookiecutter_path, "cookiecutter.json")) as fp:
    cc = json.load(fp)


def generate():
    for platform in cc["platform"]:
        for py in cc["python"]:
            cookiecutter(
                template=cookiecutter_path,
                no_input=True,
                overwrite_if_exists=True,
                output_dir=BUILD,
                extra_context={
                    "python": py,
                    "platform": platform
                }
            )


def build():
    builders = list(glob(join(BUILD, "*", "build.sh")))
    for builder in builders:
        build_dir = dirname(builder)

        if not len(list(glob(join(build_dir, "psi*.sh")))):
            p = subprocess.Popen(["bash", builder], cwd=dirname(builder))
            p.wait()

if __name__ == "__main__":
    generate()
    build()
