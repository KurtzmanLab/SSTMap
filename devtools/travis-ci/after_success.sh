#!/bin/bash
# Must be invoked with $PACKAGENAME

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [ "$TRAVIS_PULL_REQUEST" = true ]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [ "$TRAVIS_BRANCH" != "master" ]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi



