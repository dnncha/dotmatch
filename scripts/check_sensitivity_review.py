#!/usr/bin/env python3
"""Exercise the actual portable reviewer in a real browser and Node's Web Crypto.

Development-only: pytest and playwright are required. This does not claim a
native matcher run, installed-wheel validation, Safari or a desktop Tauri test.
Browser navigation can be disabled by managed environments: DOM mode still
checks layout/interactions, while Node checks the production attachment logic.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def fixtures():
    spec = importlib.util.spec_from_file_location("review_contract_tests", ROOT / "python/tests/test_sensitivity_review.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NODE_TESTS = r'''
const assert = require("node:assert/strict"), fs = require("node:fs"), vm=require("node:vm");
const {webcrypto}=require("node:crypto");
const model=JSON.parse(fs.readFileSync(process.argv[2],"utf8")), contents=fs.readFileSync(process.argv[3]);
const script=fs.readFileSync(process.argv[4],"utf8");
const nodes=new Map();const node=id=>{if(!nodes.has(id))nodes.set(id,{textContent:"",value:"",disabled:false,open:false,files:[],classList:{add(){},remove(){}},replaceChildren(){},append(){}});return nodes.get(id);};
const context={console,TextDecoder,TextEncoder,Uint8Array,ArrayBuffer,Map,Set,JSON,Number,Math,Error,Blob,File,setTimeout,clearTimeout,crypto:webcrypto,document:{getElementById:node},};context.globalThis=context;
vm.createContext(context);vm.runInContext(script,context);const api=context.__reviewTest;let count=0;
function test(name,fn){fn();count++;console.log("PASS "+name);}
function plain(value){return JSON.parse(JSON.stringify(value));}
const clone=()=>plain(model);
api.load(model);
test("embedded data reconciles",()=>api.validateEmbedded(model));
test("tampered count fails",()=>{const x=clone();x.guides[0][2]++;assert.throws(()=>api.validateEmbedded(x));});
test("unknown schema fails",()=>{const x=clone();x.schema_version="future";assert.throws(()=>api.validateEmbedded(x));});
test("duplicate guide fails",()=>{const x=clone();x.guides[1][0]=x.guides[0][0];assert.throws(()=>api.validateEmbedded(x));});
test("transition corruption fails",()=>{const x=clone();x.transitions[0][4]++;assert.throws(()=>api.validateEmbedded(x));});
test("TSV quoted fields tabs newlines CRLF",()=>assert.deepEqual(plain(api.parseTSV('a\tb\r\n"x\ty"\t"z\nq"\r\n')),[["a","b"],["x\ty","z\nq"]]));
test("TSV escaped quotes and BOM",()=>assert.deepEqual(plain(api.parseTSV('\ufeffa\tb\n"x""y"\tz\n')),[["a","b"],['x"y',"z"]]));
test("TSV unterminated quote rejected",()=>assert.throws(()=>api.parseTSV('a\tb\n"x\ty')));
test("TSV trailing garbage rejected",()=>assert.throws(()=>api.parseTSV('a\tb\n"x"oops\ty\n')));
test("TSV overlong field rejected",()=>assert.throws(()=>api.parseTSV("x".repeat(16385))));
const records=api.validateReadChanges(contents.toString("utf8"));
test("all five recorded changes reconcile",()=>assert.equal(records.length,5));
test("exact/radius ambiguity has rule explanation",()=>assert.match(api.decisionReason(records[0]),/another target/));
test("record duplicate ordinal rejected",()=>{const text=contents.toString("utf8").replace(/\n2\t/,"\n1\t");assert.throws(()=>api.validateReadChanges(text));});
test("changed row missing rejected",()=>assert.throws(()=>api.validateReadChanges(contents.toString("utf8").split("\n").slice(0,-2).join("\n")+"\n")));
test("read IDs may repeat",()=>{const text=contents.toString("utf8").replace(/synthetic_record_\d+/g,"repeat");assert.equal(api.validateReadChanges(text).length,5);});
test("unknown target rejected",()=>assert.throws(()=>api.validateReadChanges(contents.toString("utf8").replace(/guide_A/g,"missing_guide"))));
test("Hamming policy contradiction rejected",()=>{const r=plain(records[0]);r.best_k1_target_id="guide_B";assert.equal(api.validPolicyCalls(r),false);});
test("non-unique target ID rejected",()=>{const lines=contents.toString("utf8").trimEnd().split("\n"),cols=lines[1].split("\t");cols[5]="guide_A";lines[1]=cols.join("\t");assert.throws(()=>api.validateReadChanges(lines.join("\n")+"\n"));});
test("changed-read aggregate mismatch rejected",()=>{const text=contents.toString("utf8").replace(/guide_A/g,"guide_C");assert.throws(()=>api.validateReadChanges(text));});
async function attachment(name,file,pattern){node("read-file").files=[file];await api.attach();assert.match(node("attachment-status").textContent,pattern);count++;console.log("PASS "+name);}
(async()=>{
 await attachment("actual Web Crypto hash and aggregate acceptance",new File([contents],"read_changes.tsv"),/5 changed records loaded/);
 assert.equal(api.state().attached.length,5);
 api.clearReads();assert.equal(api.state().attached,null);assert.equal(node("read-file").value,"");count++;
 const changed=Buffer.from(contents);changed[changed.length-2]^=1;
 await attachment("same-size SHA mismatch refused",new File([changed],"read_changes.tsv"),/SHA-256 mismatch/);assert.equal(api.state().attached,null);
 await attachment("wrong size refused",new File([contents,Buffer.from("x")],"read_changes.tsv"),/size does not match/);
 context.crypto={};await attachment("unavailable Web Crypto fails closed",new File([contents],"read_changes.tsv"),/cannot check SHA-256/);context.crypto=webcrypto;
 let resolve;const slow={size:contents.length,arrayBuffer:()=>new Promise(r=>{resolve=r;})};node("read-file").files=[slow];const pending=api.attach();api.clearReads();resolve(contents.buffer.slice(contents.byteOffset,contents.byteOffset+contents.byteLength));await pending;assert.equal(api.state().attached,null);assert.match(node("attachment-status").textContent,/cleared/);count++;console.log("PASS clearing an in-flight attachment wins");
 node("read-file").files=[];await api.attach();assert.equal(api.state().attached,null);count++;
 console.log(JSON.stringify({node_contract_checks:count}));
})().catch(error=>{console.error(error);process.exitCode=1;});
'''


def node_checks(t, work):
    bundle = work / "node-fixture"
    summary = t.fixture(bundle)
    model = t.review.build_review_data(summary, bundle)
    (work / "model.json").write_text(json.dumps(model))
    script = t.review.JS
    boundary = 'try{\n d=JSON.parse($("review-data").textContent);'
    assert script.count(boundary) == 1
    # Use the exact production parser/validator/attachment functions, while
    # exposing them to the test harness instead of bootstrapping browser UI.
    expose = '''globalThis.__reviewTest={validateEmbedded,parseTSV,validateReadChanges,decisionReason,validPolicyCalls,attach,clearReads,
load(data){d=data;s=d.summary;},state(){return {attached};}}; return;
try{\n d=JSON.parse($("review-data").textContent);'''
    (work / "production-review.js").write_text(script.replace(boundary, expose))
    (work / "node-checks.cjs").write_text(NODE_TESTS)
    result = subprocess.run(["node", str(work / "node-checks.cjs"), str(work / "model.json"), str(bundle / "read_changes.tsv"), str(work / "production-review.js")], capture_output=True, text=True)
    print(result.stdout, end="")
    if result.returncode:
        raise RuntimeError(result.stderr)
    return json.loads(result.stdout.strip().splitlines()[-1])


def browser_checks(t, work, executable=None, screenshot_dir=None, file_mode=False, engine="chromium"):
    from playwright.sync_api import sync_playwright, expect
    bundle = work / "browser-fixture"
    summary = t.fixture(bundle, extra_zero_guides=110)
    model = t.review.build_review_data(summary, bundle)
    document = t.review.render_sensitivity_report(summary, bundle)
    tests, errors, requests = 0, [], []
    def check(condition, message):
        nonlocal tests
        assert condition, message
        tests += 1
        print("PASS " + message)
    with sync_playwright() as p:
        kwargs = {"headless": True}
        if executable: kwargs["executable_path"] = executable
        browser = getattr(p, engine).launch(**kwargs)
        context = browser.new_context(viewport={"width":1440,"height":1000}, accept_downloads=True)
        page = context.new_page()
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.on("request", lambda r: requests.append(r.url) if r.url.startswith(("http:","https:","ws:","wss:")) else None)
        if file_mode:
            report_path=bundle/"report.html";report_path.write_text(document)
            page.goto(report_path.as_uri())
        else:
            page.set_content(document)
        native_crypto=page.evaluate("!!globalThis.crypto?.subtle")
        if file_mode and not native_crypto:
            raise AssertionError("File-mode gate requires real browser Web Crypto; do not replace it with a digest double")
        expect(page.locator("#app")).to_be_visible()
        check(page.locator("#headline").inner_text()=="Assignment sensitivity", "descriptive report title")
        check(page.locator("#hero-description").inner_text()=="3 of 115 guide counts differ. Unique assignments: 3 → 3 (change: 0).", "same-total comparison reports actual guide differences")
        check(page.locator("#outcomes-title").inner_text()=="Read outcomes", "descriptive read-outcome heading")
        check(page.locator("#guide-rows tr").count()==3, "default changed-guide list")
        page.select_option("#right","best_k1")
        check(page.locator("#metric-delta").inner_text()=="+2", "policy switch updates actual counts")
        check("Unique assignments: 3 → 5 (change: +2)." in page.locator("#hero-description").inner_text(), "comparison summary updates with selected policy")
        page.click("#swap")
        check(page.locator("#metric-delta").inner_text()=="−2", "reverse comparison changes delta sign")
        check("Unique assignments: 5 → 3 (change: −2)." in page.locator("#hero-description").inner_text(), "reversed summary agrees with signed delta")
        page.click("#swap");page.select_option("#right","radius_k1")
        page.click("#reset")
        check(page.locator("#guide-rows tr").count()==50, "DOM bounded to 50 rows")
        page.click("#next")
        check("51–100 of 115" in page.locator("#page-label").inner_text(), "pagination reaches beyond original 50-row limit")
        page.fill("#search","zero_000109")
        expect(page.locator("#page-label")).to_contain_text("1–1 of 1")
        check("zero_000109" in page.locator("#guide-rows").inner_text(), "full-library search reaches last guide")
        page.fill("#search","no_such_guide")
        expect(page.locator("#guide-rows")).to_contain_text("No guides match")
        check(page.locator("#export").is_disabled(), "empty filter cannot export misleading file")
        page.click("#reset")
        with page.expect_download() as download:
            page.click("#export")
        path = work / "all-guides.tsv";download.value.save_as(path)
        check(len(path.read_text().splitlines())==116, "export contains every filtered row, not only visible page")
        with page.expect_download() as svg_download:
            page.click("#export-figure")
        figure_path = work / "figure.svg";svg_download.value.save_as(figure_path)
        import xml.etree.ElementTree as ET
        figure = ET.fromstring(figure_path.read_text())
        check(figure.tag == "{http://www.w3.org/2000/svg}svg", "SVG figure is valid self-contained XML")
        figure_text = " ".join(figure.itertext())
        check("12 of 115 matching" in figure_text and "same linear" not in figure_text, "SVG explicitly labels bounded filtered guide slice")
        check("synthetic_record_" not in figure_text, "SVG excludes optional read IDs")
        page.fill("#search","guide_A")
        expect(page.locator("#page-label")).to_contain_text("1–1 of 1")
        target = page.locator("#guide-rows .target").first
        target.click()
        expect(page.locator("#guide-dialog")).to_be_visible()
        check(page.locator("#detail-counts strong").all_inner_texts()==["1","0","1"], "guide inspector shows all three recorded counts")
        check("cannot be inferred" in page.locator("#guide-dialog").inner_text(), "absent base-level evidence is explicit")
        page.keyboard.press("Escape")
        expect(target).to_be_focused()
        check(not page.locator("#guide-dialog").is_visible(), "Escape closes inspector and restores focus")
        cell=page.get_by_role("button",name="Unique to Ambiguous: 2 reads",exact=True)
        cell.click()
        expect(page.get_by_role("button",name="Unique to Ambiguous: 2 reads",exact=True)).to_be_focused()
        check("2 reads: Unique → Ambiguous" in page.locator("#transition-selection").inner_text(), "transition selection and keyboard focus survive rerender")
        # about:blank has no secure-context Web Crypto. This MUST refuse rather
        # than silently mark an attachment verified. Node exercises success.
        page.set_input_files("#read-file",str(bundle / "read_changes.tsv"))
        if native_crypto:
            expect(page.locator("#attachment-status")).to_contain_text("5 changed records loaded")
            check("synthetic_record_1" in page.locator("#transition-reads").inner_text(), "native browser Web Crypto accepts the original attachment")
        else:
            expect(page.locator("#attachment-status")).to_contain_text("cannot check SHA-256")
            check("synthetic_record_" not in page.locator("#app").inner_text(), "unsupported browser hash context exposes no read IDs")
        page.click("#clear-reads")
        check("cleared" in page.locator("#attachment-status").inner_text(), "attachment clear state")
        # DOM acceptance/privacy checks use a DIGEST TEST DOUBLE because
        # this managed browser cannot navigate to file:// or localhost.
        # Actual Web Crypto and hash rejection are exercised separately in Node.
        if not native_crypto:
            page.expose_function("__test_digest", lambda values: list(hashlib.sha256(bytes(values)).digest()))
            page.evaluate("""() => Object.defineProperty(crypto, 'subtle', {configurable:true,
                value:{digest: async (_algorithm, bytes) => new Uint8Array(await window.__test_digest(Array.from(new Uint8Array(bytes)))).buffer}})""")
        prefix="native browser: " if native_crypto else "digest double: "
        page.set_input_files("#read-file",str(bundle / "read_changes.tsv"))
        expect(page.locator("#attachment-status")).to_contain_text("5 changed records loaded")
        check("synthetic_record_1" in page.locator("#transition-reads").inner_text(), prefix+"accepted records reach selected transition")
        page.locator("#guide-rows .target").first.click()
        check("synthetic_record_1" in page.locator("#detail-reads").inner_text(), prefix+"guide inspector receives relevant read evidence")
        page.locator("#detail-reads summary").first.click()
        check("another target" in page.locator("#detail-reads").inner_text(), prefix+"recorded-policy rationale is inspectable")
        page.keyboard.press("Escape")
        page.emulate_media(media="print")
        check(all(not table.is_visible() for table in page.locator(".read-table").all()), prefix+"print excludes attached read-ID tables")
        page.emulate_media(media="screen")
        with page.expect_download() as private_export:
            page.click("#export")
        private_path=work/"private-export.tsv";private_export.value.save_as(private_path)
        check("synthetic_record_" not in private_path.read_text(), prefix+"count export excludes attached read IDs")
        page.click("#clear-reads")
        check("synthetic_record_" not in page.locator("body").text_content(), prefix+"clear removes IDs even from closed dialog DOM")
        page.set_input_files("#read-file",str(bundle / "read_changes.tsv"))
        expect(page.locator("#attachment-status")).to_contain_text("5 changed records loaded")
        wrong=work/"wrong.tsv";raw=bytearray((bundle/"read_changes.tsv").read_bytes());raw[-2]^=1;wrong.write_bytes(raw)
        page.set_input_files("#read-file",str(wrong))
        expect(page.locator("#attachment-status")).to_contain_text("SHA-256 mismatch")
        check("synthetic_record_" not in page.locator("body").text_content(), prefix+"replacing evidence with wrong file removes stale IDs")
        page.set_input_files("#read-file",[])
        check("cleared" in page.locator("#attachment-status").inner_text(), prefix+"empty file selection clears previous evidence")
        page.click("#reset")
        if screenshot_dir:
            screenshot_dir.mkdir(parents=True, exist_ok=True)
            page.screenshot(path=str(screenshot_dir / "desktop.png"),full_page=True)
        for width in (390,768,1280):
            page.set_viewport_size({"width":width,"height":844})
            check(not page.evaluate("document.documentElement.scrollWidth > innerWidth"),f"no page overflow at {width}px")
        page.set_viewport_size({"width":390,"height":844})
        check(page.locator(".guide-table").bounding_box()["width"]<=page.locator("#guide-scroll").bounding_box()["width"], "mobile comparison fits reference, alternative and delta")
        check(page.locator(".matrix").bounding_box()["width"]<360, "all four transition columns fit mobile")
        if screenshot_dir: page.screenshot(path=str(screenshot_dir / "mobile.png"),full_page=True)
        page.emulate_media(media="print")
        check(page.locator(".print-note").is_visible(), "print discloses current-page scope")
        check(page.locator("#read-file").is_hidden(), "print excludes private attachment controls")
        page.emulate_media(media="screen")
        check(not errors, "no JavaScript exceptions")
        check(not requests, "zero outgoing requests during browser checks")
        unsafe=work/"unsafe-fixture"
        unsafe_summary=t.fixture(unsafe,identifiers=['</script><script>alert(1)</script>','b','c','d','e'])
        page.set_content(t.review.render_sensitivity_report(unsafe_summary,unsafe))
        expect(page.locator("#app")).to_be_visible()
        check("<script>alert(1)</script>" in page.locator("#guide-rows").inner_text(), "malicious-looking ID rendered as text")
        bad=model.copy();bad["schema_version"]="future"
        malformed=document.replace('"schema_version":"dotmatch.review.v1"','"schema_version":"future"',1)
        page.set_content(malformed)
        check(page.locator("#app").is_hidden() and page.locator("#boot-error").is_visible(), "unknown embedded schema refuses interactive display")
        nojs=browser.new_context(java_script_enabled=False)
        fallback=nojs.new_page();fallback.set_content(document)
        check(fallback.locator("#fallback").is_visible() and "Recorded outcomes" in fallback.inner_text("body"), "no-JavaScript static snapshot available")
        nojs.close();context.close();browser.close()
    return {"browser_engine":engine,"browser_version":browser.version,"browser_dom_checks":tests,"javascript_errors":errors,"outgoing_requests":requests,
            "digest_test_double_dom_checks":0 if native_crypto else 8,
            "browser_scope":f"Local-file Playwright {engine} with native Web Crypto; branded Safari and Tauri not covered" if file_mode else f"DOM-loaded {engine}; 8 accepted-evidence DOM checks use a Python digest test double. Actual Web Crypto exercised in Node. File/HTTP navigation, browser secure-context acceptance, Safari and Tauri not covered"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--browser", choices=("chromium", "firefox", "webkit"), default="chromium")
    parser.add_argument("--chromium",help="Optional system Chromium executable")
    parser.add_argument("--screenshots",type=Path)
    parser.add_argument("--file-mode",action="store_true",help="Require local-file navigation and real browser Web Crypto (CI acceptance gate)")
    parser.add_argument("--result",type=Path)
    args=parser.parse_args()
    if args.chromium and args.browser != "chromium":
        parser.error("--chromium is only valid with --browser chromium")
    t=fixtures()
    with tempfile.TemporaryDirectory(prefix="dotmatch-review-") as tmp:
        work=Path(tmp)
        result={**node_checks(t,work),**browser_checks(t,work,args.chromium,args.screenshots,args.file_mode,args.browser)}
    text=json.dumps(result,indent=2)
    print(text)
    if args.result: args.result.write_text(text+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
