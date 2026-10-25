//! BAM input (test path): decodes records back to SAM text lines so the rest
//! of the pipeline has a single code path.

use std::io::BufRead;

use crate::error::{Result, input};
use crate::input::sam::LineSource;

pub fn open(_r: Box<dyn BufRead + Send>) -> Result<Box<dyn LineSource>> {
    Err(input("BAM input is not implemented yet"))
}
