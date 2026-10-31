"""Focused source-only tests. No native engine, network or external test fixtures."""
from __future__ import annotations
import copy
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "python/assaycode/evidence_passport.py"
spec = importlib.util.spec_from_file_location("evidence_passport_under_test", MODULE)
assert spec and spec.loader
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


class PassportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reads = b"@read\nACGT\n+\nIIII\n"
        self.reference = b"id\tsequence\ng1\tACGT\n"
        (self.root / "reads.fastq").write_bytes(self.reads)
        (self.root / "reference.tsv").write_bytes(self.reference)
        self.project = {"schema": "dotmatch.project/v1", "project_id": "test-project",
                        "question": "Count synthetic known targets.", "study_design": "One synthetic unit, no inference.",
                        "analysis_method": "Exact matching.", "reporting_limits": "No biological conclusion.",
                        "samples": [{"id": "s1", "biological_unit": "b1", "condition": "test"}],
                        "files": [
                            {"id": "r1", "path": "reads.fastq", "role": "reads", "sample_id": "s1",
                             "expected_sha256": hashlib.sha256(self.reads).hexdigest()},
                            {"id": "ref", "path": "reference.tsv", "role": "reference", "version": "v1",
                             "expected_sha256": hashlib.sha256(self.reference).hexdigest()}]}
        self.save_project()

    def save(self, path, value):
        (self.root / path).write_text(json.dumps(value), encoding="utf-8")

    def save_project(self):
        self.save("project.json", self.project)

    def approved(self):
        report = p.inspect_project(self.root)
        contract = p.approve_contract(p.draft_contract(report), report, "Example Reviewer", "All six decisions reviewed.")
        self.save("analysis-contract.json", contract)
        return report, contract

    def passport(self):
        self.approved()
        session = p.Session(self.root)
        session.compile()
        return session, session.export()

    def test_clean_files_require_human_review(self):
        report = p.inspect_project(self.root)
        self.assertEqual(report["status"], "needs_human_review")
        self.assertEqual(report["files"][0]["size_bytes"], len(self.reads))
        self.assertFalse(report["files"][0]["qc"]["whole_fastq_validated"])

    def test_wrong_reference_refused(self):
        (self.root / "reference.tsv").write_text("wrong library")
        report = p.inspect_project(self.root)
        self.assertEqual(report["status"], "refused")
        self.assertEqual(report["findings"][0]["code"], "FILE_IDENTITY_MISMATCH")

    def test_missing_file_refused(self):
        (self.root / "reads.fastq").unlink()
        self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_malformed_fastq_refused_even_with_matching_hash(self):
        malformed = b"@read\nACGT\n+\nIII\n"
        (self.root / "reads.fastq").write_bytes(malformed)
        self.project["files"][0]["expected_sha256"] = hashlib.sha256(malformed).hexdigest()
        self.save_project()
        self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_empty_fastq_refused(self):
        (self.root / "reads.fastq").write_bytes(b"")
        self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_gzip_reads(self):
        data = gzip.compress(self.reads, mtime=0)
        (self.root / "reads.fastq.gz").write_bytes(data)
        self.project["files"][0].update(path="reads.fastq.gz", expected_sha256=hashlib.sha256(data).hexdigest())
        self.save_project()
        self.assertEqual(p.inspect_project(self.root)["status"], "needs_human_review")

    def test_corrupt_gzip_refused(self):
        (self.root / "reads.fastq.gz").write_bytes(b"broken")
        self.project["files"][0]["path"] = "reads.fastq.gz"
        self.save_project()
        self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_prefix_is_explicitly_bounded(self):
        data = self.reads * (p.MAX_RECORDS + 1) + b"malformed tail\n"
        (self.root / "reads.fastq").write_bytes(data)
        self.project["files"][0]["expected_sha256"] = hashlib.sha256(data).hexdigest()
        self.save_project()
        report = p.inspect_project(self.root)
        self.assertEqual(report["status"], "needs_human_review")
        self.assertEqual(report["files"][0]["qc"]["records_checked"], p.MAX_RECORDS)
        self.assertFalse(report["files"][0]["qc"]["whole_fastq_validated"])

    def test_size_budget_refuses_not_truncates(self):
        with patch.object(p, "MAX_FILE", 1):
            self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_prefix_budget_refuses(self):
        with patch.object(p, "MAX_PREFIX_BYTES", 1):
            self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_paths_are_confined(self):
        for path in ("../outside", "/etc/passwd", "a/../b", "./reads.fastq", "a//b", "C:\\file", "a\\b"):
            with self.subTest(path=path), self.assertRaises(p.GateError):
                p.relative_path(path)

    def test_leaf_symlink_refused(self):
        (self.root / "reads.fastq").unlink()
        (self.root / "reads.fastq").symlink_to(self.root / "reference.tsv")
        self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_parent_symlink_refused(self):
        (self.root / "linked").symlink_to(self.root, target_is_directory=True)
        self.project["files"][0]["path"] = "linked/reads.fastq"
        self.save_project()
        self.assertEqual(p.inspect_project(self.root)["status"], "refused")

    def test_duplicate_paths_and_ids_rejected(self):
        for field in ("id", "path"):
            candidate = copy.deepcopy(self.project)
            candidate["files"][1][field] = candidate["files"][0][field]
            with self.subTest(field=field), self.assertRaises(p.GateError):
                p.validate_project(candidate)

    def test_missing_reference_version_rejected(self):
        del self.project["files"][1]["version"]
        with self.assertRaises(p.GateError):
            p.validate_project(self.project)

    def test_unknown_sample_rejected(self):
        self.project["files"][0]["sample_id"] = "not-declared"
        with self.assertRaises(p.GateError):
            p.validate_project(self.project)

    def test_project_unknown_fields_rejected(self):
        self.project["allow_unsafe"] = True
        with self.assertRaises(p.GateError):
            p.validate_project(self.project)

    def test_contract_has_exactly_six_unresolved_decisions(self):
        contract = p.draft_contract(p.inspect_project(self.root))
        self.assertEqual(set(contract["decisions"]), set(p.DECISIONS))
        self.assertTrue(all(d["confirmed"] is False for d in contract["decisions"].values()))
        with self.assertRaises(p.GateError):
            p.Session(self.root).compile()

    def test_unresolved_decision_cannot_compile_even_if_rehashed(self):
        report, contract = self.approved()
        for name in p.DECISIONS:
            candidate = copy.deepcopy(contract)
            candidate["decisions"][name]["confirmed"] = False
            candidate["approval"]["contract_sha256"] = p.digest(p.contract_body(candidate))
            with self.subTest(name=name), self.assertRaises(p.GateError):
                p.compilation_plan(candidate, report)

    def test_empty_reviewer_or_rationale_rejected(self):
        report = p.inspect_project(self.root)
        for reviewer, reason in (("", "reason"), ("Reviewer", " ")):
            with self.subTest(reviewer=reviewer), self.assertRaises(p.GateError):
                p.approve_contract(p.draft_contract(report), report, reviewer, reason)

    def test_approval_cannot_override_bad_reference(self):
        (self.root / "reference.tsv").write_text("wrong")
        report = p.inspect_project(self.root)
        with self.assertRaises(p.GateError):
            p.approve_contract(p.draft_contract(report), report, "Reviewer", "override")

    def test_contract_changes_invalidate_approval(self):
        report, contract = self.approved()
        contract["decisions"]["analysis_method"]["rationale"] = "Changed later"
        with self.assertRaises(p.GateError):
            p.compilation_plan(contract, report)

    def test_changed_project_invalidates_contract(self):
        self.approved()
        self.project["study_design"] = "New biological unit declaration"
        self.save_project()
        with self.assertRaises(p.GateError):
            p.Session(self.root).compile()

    def test_changed_file_invalidates_compilation(self):
        self.approved()
        (self.root / "reads.fastq").write_bytes(self.reads + self.reads)
        with self.assertRaises(p.GateError):
            p.Session(self.root).compile()

    def test_compile_is_not_workflow_execution(self):
        self.approved()
        plan = p.Session(self.root).compile()
        self.assertEqual(plan["execution_status"], "not_executed")
        self.assertEqual(plan["engine"], "external_adapter_required")

    def test_round_trip_and_independent_pin(self):
        _, passport = self.passport()
        verdict = p.verify_passport(p.loads(json.dumps(passport)), passport["integrity"]["payload_sha256"])
        self.assertEqual(verdict["integrity"], "verified_against_independent_digest")
        self.assertFalse(verdict["identity_authenticated"])
        self.assertFalse(verdict["external_file_bytes_rechecked"])

    def test_no_pin_does_not_claim_authenticity(self):
        _, passport = self.passport()
        self.assertEqual(p.verify_passport(passport)["integrity"], "internally_consistent_only")

    def test_tamper_and_rehash_still_fails_independent_pin(self):
        _, passport = self.passport()
        changed = copy.deepcopy(passport["payload"])
        changed["agent"]["name"] = "another-client"
        with self.assertRaises(p.GateError):
            p.verify_passport(p.seal(changed), passport["integrity"]["payload_sha256"])

    def test_simple_payload_tampering_rejected(self):
        _, passport = self.passport()
        passport["payload"]["agent"]["name"] = "changed"
        with self.assertRaises(p.GateError):
            p.verify_passport(passport)

    def test_event_reordering_rejected(self):
        _, passport = self.passport()
        payload = passport["payload"]
        payload["events"].reverse()
        with self.assertRaises(p.GateError):
            p.verify_passport(p.seal(payload))

    def test_false_ready_claim_rejected(self):
        session = p.Session(self.root)
        passport = session.export()
        passport["payload"]["status"] = "ready_for_workflow_compilation"
        with self.assertRaises(p.GateError):
            p.verify_passport(p.seal(passport["payload"]))

    def test_cross_project_contract_rejected(self):
        report, contract = self.approved()
        contract["project_id"] = "other-project"
        with self.assertRaises(p.GateError):
            p.compilation_plan(contract, report)

    def test_historical_refusal_preserved_after_repair(self):
        session = p.Session(self.root)
        (self.root / "reference.tsv").write_text("wrong")
        session.inspect()
        (self.root / "reference.tsv").write_bytes(self.reference)
        self.approved()
        session.compile()
        passport = session.export()
        p.verify_passport(passport)
        self.assertTrue(any(e["kind"] == "refusal" for e in passport["payload"]["events"]))
        self.assertEqual(len(passport["payload"]["reports"]), 2)

    def test_mutation_after_success_exports_refusal_and_retains_old_approval(self):
        session, _ = self.passport()
        (self.root / "reference.tsv").write_text("wrong")
        passport = session.export()
        p.verify_passport(passport)
        self.assertEqual(passport["payload"]["status"], "refused")
        self.assertIsNone(passport["payload"]["workflow_compilation_input"])
        self.assertTrue(any(e["kind"] == "approval_observed" for e in passport["payload"]["events"]))

    def test_revoked_approval_downgrades_export(self):
        session, _ = self.passport()
        (self.root / "analysis-contract.json").unlink()
        passport = session.export()
        p.verify_passport(passport)
        self.assertEqual(passport["payload"]["status"], "needs_human_review")

    def test_fresh_named_approval_does_not_reuse_old_plan(self):
        session, _ = self.passport()
        report = p.inspect_project(self.root)
        new = p.approve_contract(p.draft_contract(report), report, "Second Reviewer", "A second review")
        self.save("analysis-contract.json", new)
        passport = session.export()
        p.verify_passport(passport)
        self.assertIsNone(passport["payload"]["workflow_compilation_input"])
        session.compile()
        p.verify_passport(session.export())

    def test_no_absolute_paths_or_raw_reads_in_export(self):
        _, passport = self.passport()
        encoded = json.dumps(passport)
        self.assertNotIn(str(self.root), encoded)
        self.assertNotIn("@read", encoded)
        self.assertNotIn("IIII", encoded)

    def test_strict_json_duplicate_keys_floats_and_nonfinite_rejected(self):
        for raw in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":1.1}', '{"a":9007199254740992}'):
            with self.subTest(raw=raw), self.assertRaises(p.GateError):
                p.loads(raw)

    def test_unicode_and_key_order_are_deterministic(self):
        self.assertEqual(p.digest({"z": [1, True, None], "a": "Donncha O’Toole"}),
                         p.digest({"a": "Donncha O’Toole", "z": [1, True, None]}))

    def test_unknown_schema_and_omitted_limits_rejected(self):
        _, passport = self.passport()
        for field, value in (("schema", "future/v9"), ("limitations", [])):
            changed = copy.deepcopy(passport["payload"])
            changed[field] = value
            with self.subTest(field=field), self.assertRaises(p.GateError):
                p.verify_passport(p.seal(changed))

    def test_cli_refusal_exit_code(self):
        (self.root / "reference.tsv").write_text("wrong")
        result = subprocess.run([sys.executable, str(MODULE), "inspect", "--root", str(self.root)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(json.loads(result.stdout)["status"], "refused")

    def test_cli_approval_requires_explicit_confirmation(self):
        self.save("analysis-contract.json", p.draft_contract(p.inspect_project(self.root)))
        result = subprocess.run([sys.executable, str(MODULE), "approve", "--root", str(self.root),
                                 "--reviewer", "Reviewer", "--rationale", "reviewed"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)

    def rpc(self, requests):
        data = "\n".join(json.dumps(item) if not isinstance(item, str) else item for item in requests) + "\n"
        result = subprocess.run([sys.executable, str(MODULE), "serve", "--root", str(self.root)],
                                input=data, capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        return [json.loads(line) for line in result.stdout.splitlines()]

    def init_requests(self):
        return [{"jsonrpc": "2.0", "id": 1, "method": "initialize",
                 "params": {"protocolVersion": "2025-11-25", "capabilities": {}, "clientInfo": {"name": "test", "version": "1"}}},
                {"jsonrpc": "2.0", "method": "notifications/initialized"}]

    def test_mcp_lifecycle_tools_and_no_approval_tool(self):
        responses = self.rpc(self.init_requests() + [{"jsonrpc": "2.0", "id": 2, "method": "tools/list"}])
        self.assertEqual(len(responses), 2)
        names = {t["name"] for t in responses[-1]["result"]["tools"]}
        self.assertEqual(names, set(p.TOOL_NAMES))
        self.assertNotIn("approve_contract", names)

    def test_mcp_preinitialization_rejected(self):
        response = self.rpc([{"jsonrpc": "2.0", "id": 1, "method": "tools/list"}])[0]
        self.assertIn("error", response)

    def test_mcp_unknown_tools_arguments_and_malformed_input(self):
        for params in ({"name": "approve_contract"}, {"name": "compile_workflow", "arguments": {"bypass": True}}):
            with self.subTest(params=params):
                response = self.rpc(self.init_requests() + [{"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": params}])[-1]
                self.assertIn("error", response)
        self.assertIn("error", self.rpc(["not json"])[0])

    def test_mcp_oversized_frame_recovery(self):
        oversized = " " * (p.MAX_JSON + 2)
        responses = self.rpc([oversized, {"jsonrpc": "2.0", "id": 2, "method": "ping"}])
        self.assertEqual(responses[0]["error"]["code"], -32700)
        self.assertEqual(responses[1]["result"], {})

    def test_mcp_refuses_unapproved_compilation(self):
        response = self.rpc(self.init_requests() + [{"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                                                   "params": {"name": "compile_workflow", "arguments": {}}}])[-1]
        self.assertTrue(response["result"]["isError"])

    def test_live_stdio_demo_end_to_end(self):
        output = self.root / "live-demo"
        result = subprocess.run([sys.executable, str(ROOT / "scripts/evidence_passport_demo.py"), "--out", str(output)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        self.assertTrue(summary["wrong_reference_refused"])
        self.assertTrue(summary["historical_refusal_retained"])
        self.assertTrue(summary["rehashed_tampering_rejected_against_pin"])
        self.assertFalse(summary["biological_workflow_executed"])


if __name__ == "__main__":
    unittest.main()
