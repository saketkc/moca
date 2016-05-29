#!/bin/bash
py.test --verbose --mpl --cov=moca --cov-config .coveragerc --cov-report term-missing
