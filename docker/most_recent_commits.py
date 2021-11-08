#!/usr/bin/env python

import json
import urllib.request


def get_json(url):
    with urllib.request.urlopen(url) as response:
        return json.loads(response.read())


branches = get_json("https://api.github.com/repos/biorack/metatlas/branches?per_page=100")
for branch_dict in branches:
    name = branch_dict["name"]
    commit_dict = get_json(f"https://api.github.com/repos/biorack/metatlas/commits/{name}")
    print(commit_dict["sha"])
