"""Regression scaffolding for the uncertainty refactor.

Two tiers (see plan):

* **Quick tier** -- small committed fixtures under
  ``pygac/tests/data/uncertainty_regression/``; runs by default during
  ``pytest`` and gates every commit of the refactor.
* **Full tier** -- opt-in end-to-end runs over the user's local
  ``$PYGAC_REGRESSION_INPUT_DIR`` of ``*LHRR*`` files; baselines stored
  under ``$PYGAC_REGRESSION_BASELINE_DIR`` (never committed).
"""
