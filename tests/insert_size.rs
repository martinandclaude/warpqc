mod common;
use common::*;

#[test]
fn insert_size_picard_rules() {
    let (r, _) = run_fixture("insert.sam");
    let i = |p: &str| format!("/insert_size/{p}");
    assert_eq!(u(&r, &i("pairs")), 5);
    close(&r, &i("median"), 130.0);
    close(&r, &i("mad"), 30.0);
    assert_eq!(u(&r, &i("mode")), 130);
    close(&r, &i("mean"), 154.0);
    close(&r, &i("sd"), 3330f64.sqrt());
    assert_eq!((u(&r, &i("p10")), u(&r, &i("p25")), u(&r, &i("p75")), u(&r, &i("p90"))), (100, 130, 160, 250));
    assert_eq!((u(&r, &i("min")), u(&r, &i("max"))), (100, 250));
    assert_eq!(u(&r, &i("rf/pairs")), 1);
    assert_eq!(u(&r, &i("tandem/pairs")), 1);
    close(&r, &i("orientation_fractions/fr"), 5.0 / 7.0);
    let h = get(&r, &i("histogram")).as_array().unwrap();
    assert_eq!(h.len(), 251);
    assert_eq!(h[130].as_u64(), Some(2));
}
