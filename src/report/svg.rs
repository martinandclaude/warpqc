//! Hand-written SVG plots for the HTML report: lines, dots and bars on one
//! linear y-axis, with short axis labels and a legend only (no titles, no
//! in-plot annotations). Colours come from CSS variables set by the page
//! (`--series-N`, `--grid`, `--axis`, `--muted`), so light and dark themes
//! swap in one place. Bars and dots carry a `<title>` tooltip, which needs
//! no JavaScript.

use std::fmt::Write;

pub const WIDTH: f64 = 640.0;
pub const HEIGHT: f64 = 300.0;
const LEFT: f64 = 60.0;
const RIGHT: f64 = 16.0;
const TOP: f64 = 28.0;
const BOTTOM: f64 = 44.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Mark {
    Line,
    Dots,
    Bars,
}

#[derive(Debug, Clone)]
pub struct Series {
    pub name: String,
    /// Categorical colour slot, 1-based.
    pub slot: usize,
    pub mark: Mark,
    /// Points in data units. For a line, `None` y values break the line.
    pub points: Vec<(f64, Option<f64>)>,
    /// Tooltip text per point (bars and dots only); empty for none.
    pub tips: Vec<String>,
    /// Hidden from the legend (e.g. the second half of a split series).
    pub legend: bool,
}

impl Series {
    pub fn new(name: &str, slot: usize, mark: Mark, points: Vec<(f64, Option<f64>)>) -> Series {
        Series {
            name: name.into(),
            slot,
            mark,
            points,
            tips: vec![],
            legend: true,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Plot {
    /// Element id, also used in the accessible label.
    pub id: String,
    pub x_label: String,
    pub y_label: String,
    pub series: Vec<Series>,
    /// Defaults to the data range (y from 0).
    pub x_range: Option<(f64, f64)>,
    pub y_range: Option<(f64, f64)>,
    /// Categorical x-axis: tick labels at x = 0, 1, 2, ...
    pub x_categories: Vec<String>,
    /// Labelled vertical separators (genome profile contig boundaries); the
    /// label is placed in the middle of the segment that starts there.
    pub x_segments: Vec<(f64, String)>,
    /// A horizontal reference line in gridline style (e.g. y = 1).
    pub y_ref: Option<f64>,
}

pub fn esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => o.push_str("&amp;"),
            '<' => o.push_str("&lt;"),
            '>' => o.push_str("&gt;"),
            '"' => o.push_str("&quot;"),
            '\'' => o.push_str("&#39;"),
            c => o.push(c),
        }
    }
    o
}

/// About `n` round ticks (1, 2 or 5 x 10^k steps) covering [lo, hi].
pub fn nice_ticks(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater)
        || !lo.is_finite()
        || !hi.is_finite()
    {
        return vec![lo];
    }
    let raw = (hi - lo) / n.max(1) as f64;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|m| m * mag)
        .find(|s| *s >= raw)
        .unwrap();
    let first = (lo / step).ceil() * step;
    let mut t = Vec::new();
    let mut v = first;
    while v <= hi + step * 1e-9 {
        t.push(if v.abs() < step * 1e-9 { 0.0 } else { v });
        v += step;
    }
    t
}

/// Compact tick label: 0, 0.5, 12, 1.5K, 2M.
pub fn fmt_tick(v: f64) -> String {
    let a = v.abs();
    let (x, suffix) = if a >= 1e6 {
        (v / 1e6, "M")
    } else if a >= 1e3 {
        (v / 1e3, "K")
    } else {
        (v, "")
    };
    let s = if (x - x.round()).abs() < 1e-9 {
        format!("{}", x.round() as i64)
    } else if (x * 10.0 - (x * 10.0).round()).abs() < 1e-9 {
        format!("{x:.1}")
    } else {
        format!("{x:.2}")
    };
    format!("{s}{suffix}")
}

fn data_range(p: &Plot) -> ((f64, f64), (f64, f64)) {
    let xs = p.series.iter().flat_map(|s| s.points.iter().map(|q| q.0));
    let ys = p
        .series
        .iter()
        .flat_map(|s| s.points.iter().filter_map(|q| q.1));
    let (mut x0, mut x1) = (f64::INFINITY, f64::NEG_INFINITY);
    for x in xs {
        x0 = x0.min(x);
        x1 = x1.max(x);
    }
    let mut y1 = f64::NEG_INFINITY;
    let mut y0 = 0f64;
    for y in ys {
        y1 = y1.max(y);
        y0 = y0.min(y);
    }
    if !p.x_categories.is_empty() {
        (x0, x1) = (-0.5, p.x_categories.len() as f64 - 0.5);
    }
    if !x0.is_finite() {
        (x0, x1) = (0.0, 1.0);
    }
    if x1 <= x0 {
        x1 = x0 + 1.0;
    }
    if !y1.is_finite() || y1 <= y0 {
        y1 = y0 + 1.0;
    }
    (
        p.x_range.unwrap_or((x0, x1)),
        p.y_range.unwrap_or((y0, y1 * 1.05)),
    )
}

pub fn render(p: &Plot) -> String {
    let ((x0, x1), (y0, y1)) = data_range(p);
    let pw = WIDTH - LEFT - RIGHT;
    let ph = HEIGHT - TOP - BOTTOM;
    let sx = |x: f64| LEFT + (x - x0) / (x1 - x0) * pw;
    let sy = |y: f64| TOP + ph - (y.clamp(y0, y1) - y0) / (y1 - y0) * ph;
    let mut o = String::new();
    let _ = write!(
        o,
        r#"<svg xmlns="http://www.w3.org/2000/svg" class="plot" id="{id}" viewBox="0 0 {WIDTH} {HEIGHT}" role="img" aria-label="{yl} by {xl}">"#,
        id = esc(&p.id),
        yl = esc(&p.y_label),
        xl = esc(&p.x_label)
    );

    // Grid and y ticks.
    for t in nice_ticks(y0, y1, 5) {
        let y = sy(t);
        let _ = write!(
            o,
            r#"<line class="grid" x1="{LEFT}" x2="{:.1}" y1="{y:.1}" y2="{y:.1}"/><text class="tick" x="{:.1}" y="{:.1}" text-anchor="end">{}</text>"#,
            LEFT + pw,
            LEFT - 6.0,
            y + 4.0,
            fmt_tick(t)
        );
    }
    if let Some(r) = p.y_ref.filter(|r| *r > y0 && *r < y1) {
        let y = sy(r);
        let _ = write!(
            o,
            r#"<line class="ref" x1="{LEFT}" x2="{:.1}" y1="{y:.1}" y2="{y:.1}"/>"#,
            LEFT + pw
        );
    }

    // X ticks: categories, segments, or numeric.
    let base = TOP + ph;
    if !p.x_categories.is_empty() {
        for (i, c) in p.x_categories.iter().enumerate() {
            let x = sx(i as f64);
            let _ = write!(
                o,
                r#"<text class="tick" x="{x:.1}" y="{:.1}" text-anchor="end" transform="rotate(-45 {x:.1} {:.1})">{}</text>"#,
                base + 14.0,
                base + 14.0,
                esc(c)
            );
        }
    } else if !p.x_segments.is_empty() {
        for (k, (start, label)) in p.x_segments.iter().enumerate() {
            let end = p.x_segments.get(k + 1).map_or(x1, |s| s.0);
            let xs = sx(*start);
            if k > 0 {
                let _ = write!(
                    o,
                    r#"<line class="grid" x1="{xs:.1}" x2="{xs:.1}" y1="{TOP}" y2="{base:.1}"/>"#
                );
            }
            // Label every segment wide enough to hold it.
            if sx(end) - xs >= 14.0 {
                let _ = write!(
                    o,
                    r#"<text class="tick" x="{:.1}" y="{:.1}" text-anchor="middle">{}</text>"#,
                    (xs + sx(end)) / 2.0,
                    base + 16.0,
                    esc(label)
                );
            }
        }
    } else {
        for t in nice_ticks(x0, x1, 8) {
            let _ = write!(
                o,
                r#"<text class="tick" x="{:.1}" y="{:.1}" text-anchor="middle">{}</text>"#,
                sx(t),
                base + 16.0,
                fmt_tick(t)
            );
        }
    }
    let _ = write!(
        o,
        r#"<line class="axis" x1="{LEFT}" x2="{:.1}" y1="{base:.1}" y2="{base:.1}"/>"#,
        LEFT + pw
    );

    // Marks.
    let bar_series = p
        .series
        .iter()
        .filter(|s| s.mark == Mark::Bars)
        .count()
        .max(1);
    let mut bar_k = 0;
    for s in &p.series {
        let class = format!("s{}", s.slot);
        match s.mark {
            Mark::Line => {
                let mut run: Vec<String> = Vec::new();
                let flush = |run: &mut Vec<String>, o: &mut String| {
                    if run.len() >= 2 {
                        let _ = write!(
                            o,
                            r#"<polyline class="line {class}" points="{}"/>"#,
                            run.join(" ")
                        );
                    }
                    run.clear();
                };
                for &(x, y) in &s.points {
                    match y {
                        Some(y) => run.push(format!("{:.1},{:.1}", sx(x), sy(y))),
                        None => flush(&mut run, &mut o),
                    }
                }
                flush(&mut run, &mut o);
            }
            Mark::Dots => {
                for (i, &(x, y)) in s.points.iter().enumerate() {
                    if let Some(y) = y {
                        let _ = write!(
                            o,
                            r#"<circle class="dot {class}" cx="{:.1}" cy="{:.1}" r="4">"#,
                            sx(x),
                            sy(y)
                        );
                        if let Some(t) = s.tips.get(i) {
                            let _ = write!(o, "<title>{}</title>", esc(t));
                        }
                        o.push_str("</circle>");
                    }
                }
            }
            Mark::Bars => {
                let slot_w = if p.x_categories.is_empty() {
                    pw / s.points.len().max(1) as f64
                } else {
                    pw / p.x_categories.len().max(1) as f64
                };
                let w = ((slot_w * 0.8 / bar_series as f64) - 2.0).clamp(1.0, 24.0);
                let offset = (bar_k as f64 - (bar_series as f64 - 1.0) / 2.0) * (w + 2.0);
                for (i, &(x, y)) in s.points.iter().enumerate() {
                    let Some(y) = y else { continue };
                    let (top, bottom) = (sy(y.max(y0)), sy(y0.max(0.0)));
                    let h = (bottom - top).max(0.0);
                    let r = (w / 2.0).min(4.0).min(h);
                    let left = sx(x) + offset - w / 2.0;
                    // Rounded data end, square at the baseline.
                    let _ = write!(
                        o,
                        r#"<path class="bar {class}" d="M{left:.1},{bottom:.1}V{:.1}Q{left:.1},{top:.1} {:.1},{top:.1}H{:.1}Q{:.1},{top:.1} {:.1},{:.1}V{bottom:.1}Z">"#,
                        top + r,
                        left + r,
                        left + w - r,
                        left + w,
                        left + w,
                        top + r,
                    );
                    if let Some(t) = s.tips.get(i) {
                        let _ = write!(o, "<title>{}</title>", esc(t));
                    }
                    o.push_str("</path>");
                }
                bar_k += 1;
            }
        }
    }

    // Axis labels.
    let _ = write!(
        o,
        r#"<text class="label" x="{:.1}" y="{:.1}" text-anchor="middle">{}</text>"#,
        LEFT + pw / 2.0,
        HEIGHT - 4.0,
        esc(&p.x_label)
    );
    let _ = write!(
        o,
        r#"<text class="label" x="14" y="{:.1}" text-anchor="middle" transform="rotate(-90 14 {:.1})">{}</text>"#,
        TOP + ph / 2.0,
        TOP + ph / 2.0,
        esc(&p.y_label)
    );

    // Legend for two or more series, top right.
    let named: Vec<&Series> = p.series.iter().filter(|s| s.legend).collect();
    let mut seen: Vec<&str> = Vec::new();
    let named: Vec<&Series> = named
        .into_iter()
        .filter(|s| {
            let new = !seen.contains(&s.name.as_str());
            seen.push(&s.name);
            new
        })
        .collect();
    if named.len() >= 2 {
        let mut x = WIDTH - RIGHT;
        let mut items = Vec::new();
        for s in named.iter().rev() {
            let w = 22.0 + 6.5 * s.name.chars().count() as f64;
            x -= w + 12.0;
            items.push((x, s));
        }
        for (x, s) in items.into_iter().rev() {
            let _ = write!(
                o,
                r#"<rect class="swatch s{}" x="{x:.1}" y="8" width="14" height="10" rx="2"/><text class="legend" x="{:.1}" y="17">{}</text>"#,
                s.slot,
                x + 20.0,
                esc(&s.name)
            );
        }
    }
    o.push_str("</svg>");
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(svg: &str) -> roxmltree::Document<'_> {
        roxmltree::Document::parse(svg).unwrap_or_else(|e| panic!("{e}\n{svg}"))
    }

    fn count(doc: &roxmltree::Document, tag: &str) -> usize {
        doc.descendants().filter(|n| n.has_tag_name(tag)).count()
    }

    #[test]
    fn line_plot_is_well_formed_and_breaks_on_none() {
        let p = Plot {
            id: "ins".into(),
            x_label: "Insert size (bp)".into(),
            y_label: "Pairs".into(),
            series: vec![Series::new(
                "FR",
                1,
                Mark::Line,
                vec![
                    (0.0, Some(1.0)),
                    (1.0, Some(3.0)),
                    (2.0, None),
                    (3.0, Some(2.0)),
                    (4.0, Some(1.0)),
                ],
            )],
            ..Default::default()
        };
        let svg = render(&p);
        let doc = parse(&svg);
        assert_eq!(count(&doc, "polyline"), 2);
        // One series: no legend.
        assert_eq!(count(&doc, "rect"), 0);
    }

    #[test]
    fn bars_with_tooltips_legend_and_escaping() {
        let mut a = Series::new(
            "Filtered",
            1,
            Mark::Bars,
            vec![(0.0, Some(30.0)), (1.0, Some(29.0))],
        );
        a.tips = vec!["chr1: 30.0X".into(), "chr2: 29.0X".into()];
        let b = Series::new(
            "Raw <&>",
            2,
            Mark::Bars,
            vec![(0.0, Some(33.0)), (1.0, Some(32.0))],
        );
        let p = Plot {
            id: "contigs".into(),
            x_label: "Contig".into(),
            y_label: "Mean depth (X)".into(),
            series: vec![a, b],
            x_categories: vec!["chr1".into(), "chr2<script>".into()],
            ..Default::default()
        };
        let svg = render(&p);
        assert!(!svg.contains("<script>"));
        let doc = parse(&svg);
        assert_eq!(count(&doc, "path"), 4);
        assert_eq!(count(&doc, "title"), 2);
        assert_eq!(count(&doc, "rect"), 2, "legend swatches");
        assert!(svg.contains("Raw &lt;&amp;&gt;"));
    }

    #[test]
    fn empty_and_single_point_plots_render() {
        for pts in [vec![], vec![(5.0, Some(2.0))]] {
            let p = Plot {
                id: "e".into(),
                series: vec![Series::new("x", 1, Mark::Dots, pts)],
                ..Default::default()
            };
            parse(&render(&p));
        }
    }

    #[test]
    fn ticks_and_labels() {
        assert_eq!(
            nice_ticks(0.0, 100.0, 5),
            vec![0.0, 20.0, 40.0, 60.0, 80.0, 100.0]
        );
        assert_eq!(
            nice_ticks(0.0, 1.0, 5),
            vec![0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0]
        );
        assert_eq!(fmt_tick(0.0), "0");
        assert_eq!(fmt_tick(0.6000000000000001), "0.6");
        assert_eq!(fmt_tick(1500.0), "1.5K");
        assert_eq!(fmt_tick(2_000_000.0), "2M");
        assert_eq!(fmt_tick(0.25), "0.25");
    }
}
