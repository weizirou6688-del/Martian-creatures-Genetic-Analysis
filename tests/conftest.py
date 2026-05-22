"""Shared pytest configuration.

Force matplotlib's non-interactive Agg backend so the chart-generation
tests run reliably in headless and CI environments (no GUI required).
"""

import matplotlib

matplotlib.use("Agg")
