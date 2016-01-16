#!/bin/bash
if [[ $TRAVIS_PULL_REQUEST == false && $TRAVIS_BRANCH == "master" ]]
then
    echo "-- pushing docs --"

    ( cd docs/_build/html
    git init
    git config user.email "travis@travis-ci.com"
    git config user.name "Travis Bot"

    ls _static
    git add .
    git commit -m "Deployed to GitHub Pages"
    echo "https://${GH_REF}"
    ghp-import -n -p -m $(GHP_MSG) -r "https://${GHTOKEN}@${GH_REF}" docs/_build/html > /dev/null 2>&1
    )
else
    echo "-- will only push docs from master --"
fi
