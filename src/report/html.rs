//! `<prefix>.qc.html`: one self-contained page rendered from the report
//! alone (so `warpqc render` can rebuild it from the JSON). No external
//! requests, no JavaScript; the full JSON is embedded for re-extraction.

use askama::Template;

use crate::report::json::Report;
use crate::report::multiqc::cumulative;
use crate::report::svg::{Mark, Plot, Series, render};

struct Row {
    label: String,
    value: String,
    status: String,
    rule: String,
}

struct PlotHtml {
    heading: String,
    svg: String,
}

#[derive(Template)]
#[template(path = "report.html")]
struct Page {
    sample: String,
    version: String,
    input: String,
    reference: String,
    reference_md5: String,
    lanes: String,
    overall_class: String,
    overall_label: String,
    summary: Vec<Row>,
    warnings: Vec<String>,
    plots: Vec<PlotHtml>,
    lane_header: Vec<&'static str>,
    lane_rows: Vec<Vec<String>>,
    definitions: Vec<(String, String)>,
    json: String,
}

fn pct(x: f64) -> String {
    format!("{:.2} %", 100.0 * x)
}

fn num(x: f64, digits: usize) -> String {
    format!("{x:.digits$}")
}

fn opt(x: Option<f64>, f: impl Fn(f64) -> String) -> String {
    x.map_or("n/a".into(), f)
}

fn thousands(n: u64) -> String {
    let s = n.to_string();
    let mut o = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            o.push(',');
        }
        o.push(c);
    }
    o
}

fn is_mito(name: &str) -> bool {
    matches!(name, "chrM" | "MT" | "chrMT" | "M")
}

/// JSON safe inside a `<script>` element: `<` becomes `<`, so no
/// `</script>` or `<!--` can appear.
pub fn script_json(r: &Report) -> String {
    serde_json::to_string(r)
        .expect("report serialises")
        .replace('<', "\\u003c")
}

/// "warn_above 0.08" shown in the unit of the value it sits next to:
/// "warn above 8.00 %" beside a percentage, "fail below 25X" beside a depth.
fn fmt_rule(rule: &str, value: &str) -> String {
    let Some((kind, t)) = rule.split_once(' ') else {
        return rule.to_string();
    };
    let Ok(t) = t.parse::<f64>() else {
        return rule.to_string();
    };
    let shown = if value.ends_with(" %") {
        pct(t)
    } else if value.ends_with('X') {
        format!("{t}X")
    } else {
        format!("{t}")
    };
    format!("{} {shown}", kind.replace('_', " "))
}

fn summary(r: &Report) -> Vec<Row> {
    let f = &r.coverage.filtered;
    let a = &r.alignment.stats;
    let d = &r.duplicates;
    let i = &r.insert_size.fr;
    let u = &r.uniformity;
    let mito = f
        .per_contig
        .iter()
        .find(|(k, _)| is_mito(k))
        .map(|(_, v)| *v);
    let mut rows: Vec<(String, &str, String)> = vec![
        (
            "coverage.filtered.mean".into(),
            "Mean depth, filtered",
            format!("{}X", num(f.mean, 2)),
        ),
        (
            "coverage.raw.mean".into(),
            "Mean depth, raw",
            opt(r.coverage.raw.as_ref().map(|x| x.mean), |v| {
                format!("{}X", num(v, 2))
            }),
        ),
        (
            "coverage.filtered.median".into(),
            "Median depth, filtered",
            format!("{}X", num(f.median, 1)),
        ),
        (
            "coverage.filtered.autosomal_mean".into(),
            "Autosomal mean depth",
            opt(f.autosomal_mean, |v| format!("{}X", num(v, 2))),
        ),
        (
            "coverage.filtered.fold80".into(),
            "Fold-80 penalty",
            opt(f.fold80, |v| num(v, 2)),
        ),
    ];
    for (t, v) in &f.pct_ge {
        rows.push((format!("coverage.filtered.pct_ge.{t}"), "", pct(*v)));
    }
    rows.sort_by_key(|(p, _, _)| {
        // Keep pct_ge rows in numeric order after the depth rows.
        p.strip_prefix("coverage.filtered.pct_ge.")
            .and_then(|t| t.parse::<u32>().ok())
            .map_or(0, |t| t + 1)
    });
    rows.extend([
        (
            "alignment.primary_mapped_fraction".into(),
            "Primary reads mapped",
            pct(a.primary_mapped_fraction),
        ),
        (
            "alignment.properly_paired_fraction".into(),
            "Properly paired",
            pct(a.properly_paired_fraction),
        ),
        (
            "alignment.mismatch_rate".into(),
            "Mismatch rate",
            opt(a.mismatch_rate, pct),
        ),
        (
            "alignment.chimeric_read_fraction".into(),
            "Mate on another contig (MAPQ ≥ 5)",
            pct(a.chimeric_read_fraction),
        ),
        (
            "alignment.off_target_fraction".into(),
            "Reads on non-analysis contigs",
            pct(a.off_target_fraction),
        ),
        ("duplicates.rate".into(), "Duplicate rate", pct(d.rate)),
        (
            "duplicates.optical_rate".into(),
            "Optical duplicate rate",
            opt(d.optical_rate, pct),
        ),
        (
            "duplicates.cross_lane_fraction".into(),
            "Duplicates spanning lanes",
            pct(d.cross_lane_fraction),
        ),
        (
            "duplicates.est_library_size".into(),
            "Estimated library size",
            d.est_library_size.map_or("n/a".into(), thousands),
        ),
        (
            "insert_size.median".into(),
            "Insert size median (FR)",
            format!("{} bp (MAD {})", num(i.median, 0), num(i.mad, 0)),
        ),
        (
            "insert_size.mean".into(),
            "Insert size mean ± SD (FR)",
            format!("{} ± {} bp", num(i.mean, 1), num(i.sd, 1)),
        ),
        ("uniformity.pmad".into(), "PMAD", opt(u.pmad, |v| num(v, 3))),
        (
            "uniformity.pmad_gc".into(),
            "PMAD, GC-corrected",
            opt(u.pmad_gc, |v| num(v, 3)),
        ),
        (
            "uniformity.evenness".into(),
            "Evenness score",
            opt(u.evenness, |v| num(v, 3)),
        ),
        (
            "uniformity.cv".into(),
            "CV of bin counts",
            opt(u.cv, |v| num(v, 3)),
        ),
    ]);
    if let Some(m) = mito {
        rows.push((
            "coverage.filtered.per_contig.chrM".into(),
            "chrM mean depth",
            format!("{}X", num(m, 1)),
        ));
    }
    // Checks on metrics not in the list above get their own rows.
    for c in &r.qc_status.checks {
        if !rows.iter().any(|(p, _, _)| *p == c.metric) {
            rows.push((c.metric.clone(), "", opt(c.value, |v| format!("{v}"))));
        }
    }
    rows.into_iter()
        .map(|(path, label, value)| {
            let check = r.qc_status.checks.iter().find(|c| c.metric == path);
            let label = if label.is_empty() {
                match path.strip_prefix("coverage.filtered.pct_ge.") {
                    Some(t) => format!("Bases at filtered depth ≥ {t}X"),
                    None => path.clone(),
                }
            } else {
                label.to_string()
            };
            let rule = check
                .and_then(|c| c.rule.as_deref())
                .map_or(String::new(), |rule| fmt_rule(rule, &value));
            Row {
                label,
                value,
                status: check.map_or(String::new(), |c| c.status.clone()),
                rule,
            }
        })
        .collect()
}

fn plots(r: &Report) -> Vec<PlotHtml> {
    let mut out = Vec::new();
    let f = &r.coverage.filtered;
    let raw = r.coverage.raw.as_ref();

    let ins: Vec<(f64, Option<f64>)> = r
        .insert_size
        .fr
        .histogram
        .iter()
        .enumerate()
        .map(|(x, &c)| (x as f64, Some(c as f64)))
        .collect();
    out.push(PlotHtml {
        heading: "Insert size".into(),
        svg: render(&Plot {
            id: "plot-insert".into(),
            x_label: "Insert size (bp)".into(),
            y_label: "Pairs".into(),
            series: vec![Series::new("FR pairs", 1, Mark::Line, ins)],
            ..Default::default()
        }),
    });

    let cum = |h: &[u64]| -> Vec<(f64, Option<f64>)> {
        cumulative(h)
            .into_iter()
            .map(|(d, p)| (d as f64, Some(p)))
            .collect()
    };
    let mut series = vec![Series::new("Filtered", 1, Mark::Line, cum(&f.histogram))];
    if let Some(rw) = raw {
        series.push(Series::new("Raw", 2, Mark::Line, cum(&rw.histogram)));
    }
    out.push(PlotHtml {
        heading: "Cumulative coverage".into(),
        svg: render(&Plot {
            id: "plot-cumulative".into(),
            x_label: "Depth (X)".into(),
            y_label: "% bases ≥ depth".into(),
            series,
            y_range: Some((0.0, 100.0)),
            ..Default::default()
        }),
    });

    // Depth histogram up to the depth that covers 99.5 % of the territory.
    let total: u64 = f.histogram.iter().sum();
    let mut acc = 0u64;
    let mut xmax = f.histogram.len().saturating_sub(1);
    for (d, &c) in f.histogram.iter().enumerate() {
        acc += c;
        if total > 0 && acc as f64 >= 0.995 * total as f64 {
            xmax = d;
            break;
        }
    }
    let xmax = xmax.max(10);
    let hist = |h: &[u64]| -> Vec<(f64, Option<f64>)> {
        let t: u64 = h.iter().sum::<u64>().max(1);
        (0..=xmax)
            .map(|d| {
                (
                    d as f64,
                    Some(100.0 * *h.get(d).unwrap_or(&0) as f64 / t as f64),
                )
            })
            .collect()
    };
    let mut series = vec![Series::new("Filtered", 1, Mark::Line, hist(&f.histogram))];
    if let Some(rw) = raw {
        series.push(Series::new("Raw", 2, Mark::Line, hist(&rw.histogram)));
    }
    out.push(PlotHtml {
        heading: "Depth distribution".into(),
        svg: render(&Plot {
            id: "plot-depth".into(),
            x_label: "Depth (X)".into(),
            y_label: "% of bases".into(),
            series,
            x_range: Some((0.0, xmax as f64)),
            ..Default::default()
        }),
    });

    // Per-contig mean depth; chrM is in the summary table (its depth would
    // flatten every other bar).
    let contigs: Vec<&String> = f.per_contig.keys().filter(|k| !is_mito(k)).collect();
    let mut ordered = contigs.clone();
    ordered.sort_by_key(|k| {
        let s = k.strip_prefix("chr").unwrap_or(k);
        match s.parse::<u32>() {
            Ok(n) => (0, n, String::new()),
            Err(_) => (1, 0, s.to_string()),
        }
    });
    let bars = |name: &str, slot: usize, m: &std::collections::BTreeMap<String, f64>| {
        let mut s = Series::new(
            name,
            slot,
            Mark::Bars,
            ordered
                .iter()
                .enumerate()
                .map(|(i, k)| (i as f64, m.get(*k).copied()))
                .collect(),
        );
        s.tips = ordered
            .iter()
            .map(|k| format!("{k} {name}: {:.2}X", m.get(*k).copied().unwrap_or(0.0)))
            .collect();
        s
    };
    let mut series = vec![bars("Filtered", 1, &f.per_contig)];
    if let Some(rw) = raw {
        series.push(bars("Raw", 2, &rw.per_contig));
    }
    out.push(PlotHtml {
        heading: "Mean depth per contig".into(),
        svg: render(&Plot {
            id: "plot-contigs".into(),
            x_label: "Contig".into(),
            y_label: "Mean depth (X)".into(),
            series,
            x_categories: ordered
                .iter()
                .map(|k| k.strip_prefix("chr").unwrap_or(k).to_string())
                .collect(),
            ..Default::default()
        }),
    });

    let u = &r.uniformity;
    let pts: Vec<(f64, Option<f64>)> = u
        .gc_curve
        .iter()
        .map(|g| (g.gc as f64, Some(g.normalized_coverage)))
        .collect();
    let mut line = Series::new("Normalised coverage", 1, Mark::Line, pts.clone());
    line.legend = false;
    let mut dots = Series::new("Normalised coverage", 1, Mark::Dots, pts);
    dots.tips = u
        .gc_curve
        .iter()
        .map(|g| {
            format!(
                "GC {} %: {:.3} ({} bins)",
                g.gc,
                g.normalized_coverage,
                thousands(g.bins)
            )
        })
        .collect();
    out.push(PlotHtml {
        heading: "GC bias".into(),
        svg: render(&Plot {
            id: "plot-gc".into(),
            x_label: "GC (%)".into(),
            y_label: "Normalised coverage".into(),
            series: vec![line, dots],
            x_range: Some((0.0, 100.0)),
            y_ref: Some(1.0),
            ..Default::default()
        }),
    });

    // Genome profile: contigs laid end to end, one line per contig.
    let mut series = Vec::new();
    let mut segments = Vec::new();
    let mut x = 0f64;
    for c in &u.profile.contigs {
        segments.push((x, c.name.strip_prefix("chr").unwrap_or(&c.name).to_string()));
        let pts: Vec<(f64, Option<f64>)> = c
            .values
            .iter()
            .enumerate()
            .map(|(i, v)| (x + i as f64, *v))
            .collect();
        let mut s = Series::new("Coverage", 1, Mark::Line, pts);
        s.legend = false;
        series.push(s);
        x += c.values.len() as f64;
    }
    out.push(PlotHtml {
        heading: format!(
            "Coverage along the genome ({} Mb windows)",
            u.profile.window_bp / 1_000_000
        ),
        svg: render(&Plot {
            id: "plot-profile".into(),
            x_label: "Autosome".into(),
            y_label: "Normalised read starts".into(),
            series,
            x_range: Some((0.0, x.max(1.0))),
            y_range: Some((0.0, 2.0)),
            x_segments: segments,
            y_ref: Some(1.0),
            ..Default::default()
        }),
    });
    out
}

pub fn render_html(r: &Report) -> String {
    let overall = r.qc_status.overall.clone();
    let lane_header = vec![
        "Lane",
        "Templates",
        "Reads",
        "Mapped",
        "Pairs",
        "Duplicates",
        "Optical",
        "Insert median",
        "Mismatch rate",
    ];
    let lane_rows = r
        .per_lane
        .iter()
        .map(|(name, l)| {
            vec![
                name.clone(),
                thousands(l.templates),
                thousands(l.reads),
                pct(l.mapped_fraction),
                thousands(l.pairs_examined),
                pct(l.duplicate_rate),
                opt(l.optical_rate, pct),
                opt(l.insert_median, |v| format!("{} bp", num(v, 0))),
                opt(l.mismatch_rate, pct),
            ]
        })
        .collect();
    let page = Page {
        sample: r.sample.clone(),
        version: r.tool.version.clone(),
        input: r.inputs.input.clone(),
        reference: r.inputs.reference.clone(),
        reference_md5: r.inputs.reference_md5.clone(),
        lanes: r.inputs.lanes.join(", "),
        overall_class: overall.clone(),
        overall_label: overall.replace('_', " "),
        summary: summary(r),
        warnings: r.warnings.clone(),
        plots: plots(r),
        lane_header,
        lane_rows,
        definitions: r
            .definitions
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect(),
        json: script_json(r),
    };
    page.render().expect("template renders")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn script_json_cannot_close_the_script_element() {
        let r = Report {
            sample: "</script><script>alert(1)</script>".into(),
            ..Default::default()
        };
        let j = script_json(&r);
        assert!(!j.contains("</"));
        let back: Report = serde_json::from_str(&j).unwrap();
        assert_eq!(back.sample, r.sample);
    }

    #[test]
    fn rules_use_the_value_unit() {
        assert_eq!(fmt_rule("warn_above 0.08", "10.02 %"), "warn above 8.00 %");
        assert_eq!(fmt_rule("fail_below 25", "6.57X"), "fail below 25X");
        assert_eq!(fmt_rule("warn_above 0.25", "0.234"), "warn above 0.25");
    }

    #[test]
    fn thousands_separators() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(1234567), "1,234,567");
        assert_eq!(thousands(999), "999");
    }
}
