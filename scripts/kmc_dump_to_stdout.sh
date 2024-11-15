#!/usr/bin/env bash

kmc_tools transform "$1" dump -s 2>/dev/null >(cat)
