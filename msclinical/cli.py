"""msclinical <tool> [args]: one command for the tools in this package."""

import importlib
import sys

TOOLS = ('edssplus', 'pirma', 'cogworse', 'neda', 'relapses', 'epochs', 'edssnoise', 'multiverse', 'tables', 'feasibility', 'table1', 'reportparse', 'sim',)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help") or argv[0] not in TOOLS:
        print("usage: msclinical <tool> [args]\ntools: " + ", ".join(TOOLS))
        return 0 if argv and argv[0] in ("-h", "--help") else 2
    return importlib.import_module(f"msclinical.{argv[0]}.cli").main(argv[1:]) or 0


if __name__ == "__main__":
    sys.exit(main())
