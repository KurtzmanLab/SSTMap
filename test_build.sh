#!/bin/bash
conda remove sstmap
conda build devtools/conda-recipe/
conda install --use-local sstmap
