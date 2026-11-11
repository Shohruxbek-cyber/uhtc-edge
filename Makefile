# Every target is defined in reproduce.py; `python reproduce.py list` shows them with their cost.
# This file only forwards, so that `make fig6` and `python reproduce.py fig6` are the same thing.
PYTHON ?= python

.PHONY: list check all
list:
	$(PYTHON) reproduce.py list

%:
	$(PYTHON) reproduce.py $@
