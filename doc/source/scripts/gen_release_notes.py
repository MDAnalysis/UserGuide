import pathlib
import os

from github import Github


def gen_release_notes(filename):
    git = Github(os.environ['GITHUB_TOKEN'])
    repo = git.get_repo('MDAnalysis/mdanalysis')
    releases = repo.get_releases()

    parent_directory = pathlib.Path(__file__).parent.parent
    parent_directory.mkdir(exist_ok=True, parents=True)
    filename = parent_directory / filename

    filetext = "# Release notes for the MDAnalysis library\n\n\n"

    # Should be ordered
    for release in repo.get_releases():
        # MDAnalysis releases always follow a tag pattern of *-release_version
        version = release.tag_name.split('-')[0]

        # Only write out version 2.x+ since those are the only ones that
        # we can guarantee similarly written notes for
        if version.split('.')[0] < 2:
            continue

        if release.body.startswith('###'):
            filetext += release.body[1:]
        else:
            filetext += release.body

        filetext += "\n\n"

    with open(filename, 'w') as f:
        f.write(filetext)


if __name__ == "__main__":
    gen_release_notes('releases.md')
