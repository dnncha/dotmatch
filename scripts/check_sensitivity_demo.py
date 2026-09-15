#!/usr/bin/env python3
"""Check the actual native-produced portable demo with real local-file Web Crypto.

Playwright WebKit on macOS is WebKit acceptance, not a branded Safari/device test.
No digest doubles or network-backed viewer assets are permitted in this gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from playwright.sync_api import expect, sync_playwright


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example", required=True, type=Path)
    parser.add_argument("--browser", choices=("chromium", "firefox", "webkit"), default="chromium")
    parser.add_argument("--result", required=True, type=Path)
    args = parser.parse_args(argv)
    example = args.example.resolve()
    errors, requests, checks = [], [], []
    with sync_playwright() as p:
        browser = getattr(p, args.browser).launch()
        context = browser.new_context(viewport={"width": 1440, "height": 1000}, accept_downloads=True)
        page = context.new_page()
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.on("request", lambda r: requests.append(r.url) if r.url.startswith(("http:", "https:", "ws:", "wss:")) else None)
        page.goto((example / "report.html").as_uri())
        expect(page.get_by_role("note", name="Synthetic demonstration")).to_be_visible()
        expect(page.locator("#headline")).to_contain_text("Same assigned total")
        assert page.evaluate("!!globalThis.crypto?.subtle")
        checks.append("labelled_native_report_with_real_browser_crypto")
        assert page.locator("#guide-rows tr").count() == 3
        page.get_by_role("button", name="guide_A", exact=True).click()
        assert page.locator("#detail-counts strong").all_inner_texts() == ["1", "0", "1"]
        page.keyboard.press("Escape")
        checks.append("native_guide_counts")
        page.set_input_files("#read-file", str(example / "bundle/read_changes.tsv"))
        expect(page.locator("#attachment-status")).to_contain_text("5 changed records loaded")
        page.get_by_role("button", name="guide_A", exact=True).click()
        expect(page.locator("#detail-reads")).to_contain_text("exact_near_a")
        page.keyboard.press("Escape")
        checks.append("native_producer_attachment_and_read_inspection")
        page.click("#clear-reads")
        assert "exact_near_a" not in page.locator("body").text_content()
        checks.append("clear_removes_native_read_identifiers")
        with page.expect_download() as event:
            page.click("#export")
        content = Path(event.value.path()).read_text()
        assert "guide_A" in content and "exact_near_a" not in content
        checks.append("filtered_export_without_read_identifiers")
        page.screenshot(path=str(args.result.with_suffix(".png")), full_page=True)
        assert not errors, errors
        assert not requests, requests
        result = {"browser_engine": args.browser, "browser_version": browser.version,
                  "checks": checks, "javascript_errors": errors, "outgoing_requests": requests,
                  "digest_test_doubles": 0, "scope": "Native-generated synthetic example; not biological validation or branded Safari acceptance"}
        context.close(); browser.close()
    args.result.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
