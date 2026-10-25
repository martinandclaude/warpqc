mod common;
use common::*;
use warpqc::report::json::Report;

fn write(sample: &str, fx: &str) -> (Report, String, tempfile::TempDir) {
    let mut c = config(fx);
    c.sample = sample.into();
    c.bin_size = 20;
    c.thresholds = Some(root().join("tests/fixtures/thresholds/qc.yaml"));
    let out = warpqc::pipeline::run(&c).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let files = warpqc::report::write_all(&out.report, &dir.path().join("S")).unwrap();
    let names: Vec<String> = files
        .iter()
        .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    assert_eq!(
        names,
        [
            "S.qc.json",
            "S.qc.html",
            "S.general_mqc.json",
            "S.insert_mqc.json",
            "S.coverage_mqc.json"
        ]
    );
    let html = std::fs::read_to_string(dir.path().join("S.qc.html")).unwrap();
    (out.report, html, dir)
}

fn between<'a>(s: &'a str, start: &str, end: &str) -> &'a str {
    let a = s.find(start).unwrap() + start.len();
    let b = a + s[a..].find(end).unwrap();
    &s[a..b]
}

#[test]
fn html_has_six_wellformed_plots_and_the_full_json() {
    let (report, html, _d) = write("S01", "coverage.sam");
    let mut ids = vec![];
    let mut rest = html.as_str();
    while let Some(i) = rest.find("<svg") {
        let end = rest[i..].find("</svg>").unwrap() + i + 6;
        let doc = roxmltree::Document::parse(&rest[i..end]).unwrap_or_else(|e| panic!("{e}"));
        ids.push(doc.root_element().attribute("id").unwrap().to_string());
        rest = &rest[end..];
    }
    assert_eq!(
        ids,
        [
            "plot-insert",
            "plot-cumulative",
            "plot-depth",
            "plot-contigs",
            "plot-gc",
            "plot-profile"
        ]
    );
    let json = between(
        &html,
        r#"<script type="application/json" id="warpqc-report">"#,
        "</script>",
    );
    assert_eq!(Report::from_json(json.as_bytes()).unwrap(), report);
    // Status from the thresholds is shown.
    assert!(html.contains(r#"<span class="badge fail">fail</span>"#));
}

#[test]
fn html_makes_no_external_requests() {
    let (_, html, _d) = write("S01", "insert.sam");
    let lower = html.to_lowercase();
    for bad in ["src=", "href=", "url(", "@import", "<link", "<script src"] {
        assert!(!lower.contains(bad), "found {bad}");
    }
    // The only URL is the SVG namespace.
    assert_eq!(
        lower.matches("http").count(),
        lower.matches("http://www.w3.org/2000/svg").count()
    );
}

#[test]
fn sample_name_is_escaped() {
    let (_, html, _d) = write("S<script>alert(1)</script>", "insert.sam");
    assert!(!html.contains("<script>alert"));
    assert!(html.contains("S&#60;script&#62;alert(1)") || html.contains("S&lt;script&gt;alert(1)"));
    assert_eq!(
        html.matches("<script").count(),
        1,
        "only the JSON script element"
    );
}

#[test]
fn per_lane_insert_median_and_mismatch_rate() {
    let (r, _) = run_fixture("insert.sam");
    close(&r, "/per_lane/FCA:1/insert_median", 130.0);
    close(&r, "/per_lane/FCA:1/mismatch_rate", 0.0);
    let (r, _) = run_fixture("flags_mixed.sam");
    // Same reads as alignment.mismatch_rate: 2 mismatches over 109 bases.
    close(&r, "/per_lane/FCA:1/mismatch_rate", 2.0 / 109.0);
}
