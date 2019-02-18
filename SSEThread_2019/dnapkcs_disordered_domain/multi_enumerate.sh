#!/bin/bash



../../../../BUILD/setup_environment.sh python modeling_enumerate_multi.py 0 100 &
../../../../BUILD/setup_environment.sh python modeling_enumerate_multi.py 1 100 &
../../../../BUILD/setup_environment.sh python modeling_enumerate_multi.py 2 100 &
../../../../BUILD/setup_environment.sh python modeling_enumerate_multi.py 3 100 &

