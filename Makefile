OCTAVE ?= octave
MATLAB ?= matlab
TEST_CODE=success = doctest({'channels', 'info', 'linalg', 'misc', 'sdp', 'states'}); exit(~success);

.PHONY: help test matlab_test

help:
	@echo Available rules:
	@echo "  test               run tests with Octave"
	@echo "  matlab_test        run tests with Matlab"
	@echo "  travis_test        run tests with Matlab"

test:
	${OCTAVE} --eval "install_quantbox; ${TEST_CODE}"

matlab_test:
	${MATLAB} -nojvm -nodisplay -nosplash -r "addpath('misc'); ${TEST_CODE}"

travis_test:
	${OCTAVE} --eval "install_quantbox; ${TEST_CODE}"
