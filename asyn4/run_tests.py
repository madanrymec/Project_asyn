"""
Zero-dependency test runner: discovers and runs every test_*() function in
tests/test_asyn4.py without needing pytest installed. Prints PASS/FAIL per
test and a summary at the end.

Run from the project root:
    python3 run_tests.py

If you have (or install) pytest, `pytest tests/` gives nicer output
(diffs on failure, -k filtering, etc) - this script is just a fallback
that needs nothing beyond numpy.
"""
import sys
import traceback
import types

# Minimal pytest.approx() shim so test_asyn4.py doesn't need pytest installed.
class _Approx:
    def __init__(self, expected, rel=1e-6, abs_=1e-12):
        self.expected = expected
        self.rel = rel
        self.abs = abs_

    def __eq__(self, other):
        try:
            iter(self.expected)
        except TypeError:
            return abs(other - self.expected) <= max(self.abs, self.rel * abs(self.expected))
        return all(
            abs(o - e) <= max(self.abs, self.rel * abs(e))
            for o, e in zip(other, self.expected)
        )

    def __radd__(self, other):
        return self

    def __repr__(self):
        return f"approx({self.expected!r})"


def _approx(expected, rel=1e-6, abs=1e-12):
    return _Approx(expected, rel=rel, abs_=abs)


class _Raises:
    def __init__(self, exc_type):
        self.exc_type = exc_type

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            raise AssertionError(f"Expected {self.exc_type.__name__} to be raised")
        if not issubclass(exc_type, self.exc_type):
            return False
        return True


def _raises(exc_type):
    return _Raises(exc_type)


_fake_pytest = types.ModuleType("pytest")
_fake_pytest.approx = _approx
_fake_pytest.raises = _raises
sys.modules.setdefault("pytest", _fake_pytest)


def main():
    import importlib
    import inspect
    import tempfile
    import pathlib

    test_mod = importlib.import_module("tests.test_asyn4")

    tests = [(name, fn) for name, fn in vars(test_mod).items()
             if name.startswith("test_") and callable(fn)]

    passed, failed = 0, 0
    for name, fn in tests:
        try:
            params = inspect.signature(fn).parameters
            if "tmp_path" in params:
                with tempfile.TemporaryDirectory() as tmpdir:
                    fn(tmp_path=pathlib.Path(tmpdir))
            else:
                fn()
            print(f"PASS  {name}")
            passed += 1
        except Exception:
            print(f"FAIL  {name}")
            traceback.print_exc()
            failed += 1

    print(f"\n{passed} passed, {failed} failed, {passed + failed} total")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
