//! `--thresholds qc.yaml`: warn/fail bounds per JSON path, evaluated into
//! `qc_status`.
//!
//! ```yaml
//! coverage.filtered.mean: {fail_below: 25, warn_below: 30}
//! duplicates.rate: {warn_above: 0.10, fail_above: 0.20}
//! ```
//!
//! A path is dot-separated JSON object keys. Paths are checked against a
//! skeleton of the report before any input is read, so a typo fails the
//! run at once (exit 2) instead of after hours of streaming. A metric that
//! comes out null is reported as "missing", which counts as a warning.
//! Failing QC never changes the exit code.

use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::error::{IoContext, Result, input};

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bounds {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fail_below: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub warn_below: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub warn_above: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fail_above: Option<f64>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Thresholds {
    /// In file order.
    pub rules: Vec<(String, Bounds)>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Check {
    pub metric: String,
    pub value: Option<f64>,
    /// "pass", "warn", "fail" or "missing".
    pub status: String,
    /// The bound that was crossed, e.g. "fail_below 25".
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub rule: Option<String>,
    pub bounds: Bounds,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct QcStatus {
    /// "pass", "warn", "fail", or "not_evaluated" without thresholds.
    pub overall: String,
    pub checks: Vec<Check>,
}

impl QcStatus {
    pub fn not_evaluated() -> QcStatus {
        QcStatus {
            overall: "not_evaluated".into(),
            checks: vec![],
        }
    }
}

fn lookup<'a>(v: &'a Value, path: &str) -> std::result::Result<&'a Value, String> {
    let mut cur = v;
    let mut walked = String::new();
    for key in path.split('.') {
        match cur {
            Value::Object(m) => match m.get(key) {
                Some(next) => cur = next,
                None => return Err(format!("no key \"{key}\" after \"{walked}\"")),
            },
            Value::Null => return Err(format!("\"{walked}\" is null in this run")),
            _ => return Err(format!("\"{walked}\" is not an object")),
        }
        if !walked.is_empty() {
            walked.push('.');
        }
        walked.push_str(key);
    }
    Ok(cur)
}

impl Thresholds {
    pub fn load(path: &Path) -> Result<Thresholds> {
        let text = std::fs::read_to_string(path).ctx(format!("reading {}", path.display()))?;
        let map: serde_norway::Mapping = serde_norway::from_str(&text)
            .map_err(|e| input(format!("{}: not a YAML mapping: {e}", path.display())))?;
        let mut rules = Vec::new();
        for (k, v) in map {
            let key = k
                .as_str()
                .ok_or_else(|| input(format!("{}: keys must be JSON paths", path.display())))?
                .to_string();
            let b: Bounds = serde_norway::from_value(v)
                .map_err(|e| input(format!("{}: {key}: {e}", path.display())))?;
            if b == Bounds::default() {
                return Err(input(format!("{}: {key} has no bounds", path.display())));
            }
            rules.push((key, b));
        }
        Ok(Thresholds { rules })
    }

    /// Checks every path against a report skeleton. Paths under `per_lane`
    /// are only known after the run and are checked then.
    pub fn validate(&self, skeleton: &Value) -> Result<()> {
        for (path, _) in &self.rules {
            if path.starts_with("per_lane.") {
                continue;
            }
            match lookup(skeleton, path) {
                Ok(Value::Number(_) | Value::Null) => {}
                Ok(_) => return Err(input(format!("threshold {path}: not a numeric metric"))),
                Err(why) => return Err(input(format!("threshold {path}: unknown metric ({why})"))),
            }
        }
        Ok(())
    }

    pub fn evaluate(&self, report: &Value) -> QcStatus {
        let mut checks = Vec::new();
        for (path, b) in &self.rules {
            let value = lookup(report, path).ok().and_then(Value::as_f64);
            let (status, rule) = match value {
                None => ("missing", None),
                Some(v) => {
                    let crossed = |bound: Option<f64>, below: bool| {
                        bound.filter(|&t| if below { v < t } else { v > t })
                    };
                    if let Some(t) = crossed(b.fail_below, true) {
                        ("fail", Some(format!("fail_below {t}")))
                    } else if let Some(t) = crossed(b.fail_above, false) {
                        ("fail", Some(format!("fail_above {t}")))
                    } else if let Some(t) = crossed(b.warn_below, true) {
                        ("warn", Some(format!("warn_below {t}")))
                    } else if let Some(t) = crossed(b.warn_above, false) {
                        ("warn", Some(format!("warn_above {t}")))
                    } else {
                        ("pass", None)
                    }
                }
            };
            checks.push(Check {
                metric: path.clone(),
                value,
                status: status.into(),
                rule,
                bounds: b.clone(),
            });
        }
        let rank = |s: &str| match s {
            "fail" => 3,
            "warn" | "missing" => 2,
            _ => 1,
        };
        let worst = checks.iter().map(|c| rank(&c.status)).max().unwrap_or(0);
        QcStatus {
            overall: match worst {
                0 => "not_evaluated",
                1 => "pass",
                2 => "warn",
                _ => "fail",
            }
            .into(),
            checks,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/thresholds")
            .join(name)
    }

    #[test]
    fn evaluates_each_bound_and_worst_wins() {
        let t = Thresholds::load(&fixture("qc.yaml")).unwrap();
        assert_eq!(t.rules.len(), 4);
        assert_eq!(t.rules[0].0, "coverage.filtered.mean");
        let report = json!({
            "coverage": {"filtered": {"mean": 0.132, "fold80": null}, "raw": {"pct_ge": {"1": 0.16}}},
            "duplicates": {"rate": 0.1},
        });
        let q = t.evaluate(&report);
        let st: Vec<(&str, &str)> = q
            .checks
            .iter()
            .map(|c| (c.metric.as_str(), c.status.as_str()))
            .collect();
        assert_eq!(
            st,
            [
                ("coverage.filtered.mean", "fail"),
                ("coverage.raw.pct_ge.1", "warn"),
                ("duplicates.rate", "pass"),
                ("coverage.filtered.fold80", "missing"),
            ]
        );
        assert_eq!(q.checks[0].rule.as_deref(), Some("fail_below 0.2"));
        assert_eq!(q.overall, "fail");
    }

    #[test]
    fn boundaries_are_inclusive_passes() {
        let t = Thresholds {
            rules: vec![(
                "x".into(),
                Bounds {
                    fail_below: Some(1.0),
                    fail_above: Some(2.0),
                    ..Default::default()
                },
            )],
        };
        assert_eq!(t.evaluate(&json!({"x": 1.0})).overall, "pass");
        assert_eq!(t.evaluate(&json!({"x": 2.0})).overall, "pass");
        assert_eq!(t.evaluate(&json!({"x": 2.5})).overall, "fail");
        assert_eq!(
            Thresholds::default().evaluate(&json!({})).overall,
            "not_evaluated"
        );
    }

    #[test]
    fn rejects_bad_files_and_paths() {
        let skeleton =
            json!({"coverage": {"filtered": {"mean": 0.0}}, "duplicates": {"rate": 0.0}});
        let typo = Thresholds::load(&fixture("typo.yaml")).unwrap();
        let e = typo.validate(&skeleton).unwrap_err();
        assert_eq!(e.exit_code(), 2);
        assert!(e.to_string().contains("meen"), "{e}");
        assert_eq!(
            Thresholds::load(&fixture("bad_key.yaml"))
                .unwrap_err()
                .exit_code(),
            2
        );
        let obj = Thresholds::load(&fixture("object.yaml")).unwrap();
        assert!(
            obj.validate(&skeleton)
                .unwrap_err()
                .to_string()
                .contains("not a numeric")
        );
        // Raw depth disabled: its paths pass through null.
        let raw = Thresholds {
            rules: vec![(
                "coverage.raw.mean".into(),
                Bounds {
                    warn_below: Some(1.0),
                    ..Default::default()
                },
            )],
        };
        assert!(raw.validate(&json!({"coverage": {"raw": null}})).is_err());
    }
}
