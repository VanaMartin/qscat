"""No-SCF controls for the fixed-QZ diffuse diagnostic and owner entrypoint."""

import ast
import hashlib
import json
import unittest
from pathlib import Path

from projects.ukrmol_co.oxygen_affinity_diffuse import INPUTS, authenticate, basis_control
from projects.ukrmol_co.oxygen_affinity_owner import container_program


class DiffuseControls(unittest.TestCase):
    def test_atomic_calculation_bytes_unchanged_from_original_qz(self):
        text = Path("projects/ukrmol_co/atomic_references.py").read_text()
        node = next(
            n
            for n in ast.parse(text).body
            if isinstance(n, ast.FunctionDef) and n.name == "calculate"
        )
        raw = ("\n".join(text.splitlines()[node.lineno - 1 : node.end_lineno]) + "\n").encode()
        self.assertEqual(
            hashlib.sha256(raw).hexdigest(),
            "a07c02721c493ff6f289c6972fb411a0694f4850ed148ec317d7e9e77a146cd0",
        )

    def test_expanded_scope_rejected(self):
        base = dict(
            source_sha256={},
            input_sha256={},
            execution_authorized=True,
            states=["O", "O-"],
            maximum_new_scf_calls=2,
            automatic_successor=False,
        )
        for change in [
            dict(states=["O"]),
            dict(maximum_new_scf_calls=4),
            dict(automatic_successor=True),
            dict(execution_authorized=False),
        ]:
            with self.assertRaises(ValueError):
                authenticate(Path("."), Path("."), base | change)

    def test_exact_parent_and_diffuse_count(self):
        inputs = json.loads(Path(INPUTS).read_text())
        self.assertEqual(
            inputs["second_diffuse_shells"],
            [
                [0, [0.0234, 1.0]],
                [1, [0.0163, 1.0]],
                [2, [0.0534, 1.0]],
                [3, [0.122, 1.0]],
                [4, [0.276, 1.0]],
            ],
        )
        try:
            import pyscf  # noqa: F401
        except ImportError:
            self.skipTest("Pinned-image integral-only AO control required before release")
        basis_control(inputs)

    def test_owner_program_compiles_without_shell_interpolation(self):
        for contract in [
            {},
            {
                "runtime_test_module": "projects.ukrmol_co.test_oxygen_affinity_diffuse",
                "runtime_science_module": "projects.ukrmol_co.oxygen_affinity_diffuse",
            },
        ]:
            program = container_program(contract)
            compile(program, "<container>", "exec")
            self.assertNotIn("shell=True", program)


if __name__ == "__main__":
    unittest.main()
