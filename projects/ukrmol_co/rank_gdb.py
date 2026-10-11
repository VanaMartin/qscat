"""Executed by GDB's Python interpreter to retain native exits or first-fault stacks."""

import json
import os
import signal
from pathlib import Path

import gdb

DIRECTORY = Path(os.environ["CO_DEBUG_DIRECTORY"])
RANK = int(os.environ.get("OMPI_COMM_WORLD_RANK", "0"))
OUTPUT = DIRECTORY / f"rank-{RANK}.json"


def save(value):
    temporary = OUTPUT.with_suffix(".tmp")
    temporary.write_text(json.dumps({"rank": RANK, **value}, indent=2) + "\n")
    temporary.replace(OUTPUT)


def exited(event):
    if OUTPUT.exists() and json.loads(OUTPUT.read_text()).get("observed_stop_signal"):
        return
    save({"inferior_exit_code": getattr(event, "exit_code", None), "observed_stop_signal": None})


def stopped(event):
    if isinstance(event, gdb.SignalEvent):
        number = int(getattr(signal, event.stop_signal))
        save(
            {
                "inferior_exit_code": None,
                "observed_stop_signal": event.stop_signal,
                "observed_signal_number": number,
                "debugger_failure_status": 128 + number,
                "scope": "Intercepted native signal; debugger status is not a native exit code",
            }
        )
        gdb.execute("thread apply all bt full")
        # The inferior is stopped at its first signal; terminate this rank so MPI aborts peers.
        gdb.execute(f"quit {128 + number}")


gdb.execute("set pagination off")
gdb.execute("set confirm off")
gdb.execute("set print elements 32")
gdb.execute("set print repeats 8")
for name in ("SIGSEGV", "SIGBUS", "SIGABRT"):
    gdb.execute(f"handle {name} stop print nopass")
gdb.events.exited.connect(exited)
gdb.events.stop.connect(stopped)
gdb.execute("run")
