#!/usr/bin/env bash

kmc_tools transform $1 dump 2>/dev/null >(cat)
