#!/usr/bin/env python3
"""Convenience entry point for verifying all theorems in this paper-companion repo.

Usage:
    python run_checks.py            # run all theorems
    python run_checks.py --verbose  # show detailed per-theorem output
    python run_checks.py --check T_Born  # run a single theorem by name

This is the standalone-repo verifier for Ledgers Entropy Time Cost (APF Paper 3).
The full APF codebase v6.8 (frozen 2026-04-18) verifies 348 checks across
335 bank-registered theorems plus apf/standalone/ — this paper's repo
verifies the 13-theorem subset directly cited by the manuscript.
"""

import sys
import argparse
import time
from apf.bank import REGISTRY, run_all, get_check


def main():
    parser = argparse.ArgumentParser(description="Verify all theorems in this paper-companion repo.")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show detailed per-theorem output")
    parser.add_argument("--check", "-c", type=str, default=None,
                        help="Run a single theorem by name (e.g., --check T_Born)")
    args = parser.parse_args()

    if args.check:
        try:
            check_fn = get_check(args.check)
        except KeyError:
            print(f"  ERROR: theorem '{args.check}' not in this repo's bank.")
            print(f"  Available: {sorted(REGISTRY.keys())}")
            sys.exit(1)
        result = check_fn()
        print(f"  {args.check}: {'PASS' if result.get('passed', True) else 'FAIL'}")
        if "key_result" in result:
            print(f"    Result: {result['key_result']}")
        if args.verbose:
            print(f"    Full result: {result}")
        sys.exit(0 if result.get('passed', True) else 1)

    start = time.time()
    results = run_all(verbose=args.verbose)
    elapsed = time.time() - start

    passed = sum(1 for r in results if r.get('passed', True))
    failed = len(results) - passed
    total = len(results)

    print()
    print("=" * 72)
    print(f"  Paper 3 ({{Ledgers Entropy Time Cost}}): {passed} passed, {failed} failed, "
          f"{total} total — verified in {elapsed:.2f}s")
    print("=" * 72)

    if failed > 0:
        print("\nFailed checks:")
        for r in results:
            if not r.get('passed', True):
                print(f"  - {r.get('name', '?')}: {r.get('error', 'no detail')}")
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
