#!/bin/bash

if [[ "${ENV}" == "DEV" ]]
then
    if [[ -f dev_requirements.txt ]]
    then
        pip install -r dev_requirements.txt
    elif [[ -f requirements-dev.txt ]]
    then
        pip install -r requirements-dev.txt
    fi
fi

exec "${@}"
