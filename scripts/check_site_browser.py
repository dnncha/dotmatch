#!/usr/bin/env python3
"""Exercise the actual static export and save reproducible viewport screenshots."""
from __future__ import annotations

import functools
import http.server
import json
import tempfile
import threading
from pathlib import Path

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]


def main():
    output = ROOT / "browser-review"
    output.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        (root / "dotmatch").symlink_to(ROOT / "out", target_is_directory=True)
        handler = functools.partial(
            http.server.SimpleHTTPRequestHandler, directory=root
        )
        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        origin = f"http://127.0.0.1:{server.server_port}"
        checks = []
        try:
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch()
                context = browser.new_context()
                page = context.new_page()
                errors = []
                page.on("pageerror", lambda error: errors.append(str(error)))
                for width in (390, 768, 1440):
                    page.set_viewport_size({"width": width, "height": 960})
                    for route in (
                        "",
                        "crispr-guide-counting/",
                        "assignment-sensitivity/",
                        "tools/library-safety/",
                    ):
                        response = page.goto(
                            origin + "/dotmatch/" + route, wait_until="networkidle"
                        )
                        assert response and response.status == 200
                        assert page.locator("h1").count() == 1
                        assert not page.evaluate(
                            "document.documentElement.scrollWidth > window.innerWidth"
                        ), (width, route, "page overflow")
                        page.screenshot(
                            path=str(
                                output
                                / f"{route.replace('/', '-') or 'home-'}{width}.png"
                            ),
                            full_page=True,
                        )
                        checks.append(
                            {"route": route, "width": width, "overflow": False}
                        )
                page.set_viewport_size({"width": 390, "height": 844})
                page.goto(origin + "/dotmatch/", wait_until="networkidle")
                for label, expected in [
                    ("Exact", "unique"),
                    ("Ambiguous", "ambiguous"),
                    ("Unmatched", "unmatched"),
                    ("Short read", "invalid"),
                ]:
                    button = page.get_by_role("button", name=label, exact=True)
                    button.click()
                    assert button.get_attribute("aria-pressed") == "true"
                    assert page.locator("[data-status]").inner_text() == expected
                menu = page.locator(".mobile-nav summary")
                menu.click()
                assert page.locator(".mobile-nav").get_attribute("open") is not None
                page.locator(".mobile-nav").get_by_role(
                    "link", name="Library checker"
                ).click()
                page.wait_for_url("**/tools/library-safety/")
                page.wait_for_load_state("networkidle")
                requests = []
                page.on(
                    "request",
                    lambda request: requests.append(
                        {
                            "url": request.url,
                            "method": request.method,
                            "body": request.post_data,
                        }
                    ),
                )
                page.get_by_role("button", name="Load synthetic example").click()
                page.get_by_role("button", name="Check complete library").click()
                page.get_by_role(
                    "heading", name="Library geometry, not an error-rate estimate."
                ).wait_for()
                assert (
                    "Complete: all 3 targets checked."
                    in page.get_by_role("status").inner_text()
                )
                with page.expect_download() as event:
                    page.get_by_role("button", name="Save complete JSON report").click()
                report = json.loads(Path(event.value.path()).read_text())
                assert report["target_count"] == 3
                assert report["ambiguous_observations"] == 2
                page.get_by_label("Target sequences", exact=True).fill(
                    "id\tsequence\na\tACGTNCGT"
                )
                assert (
                    page.get_by_role(
                        "heading", name="Library geometry, not an error-rate estimate."
                    ).count()
                    == 0
                )
                page.get_by_role("button", name="Check complete library").click()
                expect(page.locator("#library-error")).to_contain_text("A, C, G and T")
                page.get_by_role("button", name="Clear", exact=True).click()
                assert (
                    page.get_by_label("Target sequences", exact=True).input_value()
                    == ""
                )
                assert (
                    requests == []
                ), "Library interactions must not initiate network requests"
                page.goto(origin + "/dotmatch/", wait_until="networkidle")
                page.keyboard.press("Tab")
                assert (
                    page.evaluate("document.activeElement.textContent")
                    == "Skip to main content"
                )
                assert not errors, errors
                checks.append(
                    {
                        "native_example_controls": "passed",
                        "mobile_navigation": "passed",
                        "local_library_export": "passed",
                        "input_error_and_stale_result": "passed",
                        "sequence_requests": requests,
                        "keyboard_skip_link": "passed",
                        "page_errors": errors,
                    }
                )
                browser.close()
        finally:
            server.shutdown()
            server.server_close()
        (output / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
        print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
