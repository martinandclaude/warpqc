//! BAM input (test path): decodes records back to SAM text lines so the rest
//! of the pipeline has a single code path. Throughput is not a goal here.

use std::io::{BufRead, Read};

use noodles::bam;
use noodles::sam::{self, alignment::io::Write as _};

use crate::error::{IoContext, Result, input};
use crate::input::sam::LineSource;

struct BamSource<R: Read> {
    reader: bam::io::Reader<noodles::bgzf::io::Reader<R>>,
    header: sam::Header,
    text: String,
    record: bam::Record,
    writer: sam::io::Writer<Vec<u8>>,
}

impl<R: Read + Send> LineSource for BamSource<R> {
    fn header(&self) -> &str {
        &self.text
    }

    fn next_line(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        if self
            .reader
            .read_record(&mut self.record)
            .ctx("reading BAM record")?
            == 0
        {
            return Ok(false);
        }
        self.writer.get_mut().clear();
        self.writer
            .write_alignment_record(&self.header, &self.record)
            .ctx("decoding BAM record")?;
        buf.extend_from_slice(self.writer.get_ref());
        Ok(true)
    }
}

pub fn open(r: Box<dyn BufRead + Send>) -> Result<Box<dyn LineSource>> {
    let mut reader = bam::io::Reader::new(r);
    let header = reader.read_header().map_err(|e| {
        input(format!(
            "input looks like BGZF but is not a readable BAM: {e}"
        ))
    })?;
    let mut w = sam::io::Writer::new(Vec::new());
    w.write_header(&header).ctx("formatting BAM header")?;
    let text = String::from_utf8_lossy(w.get_ref()).into_owned();
    Ok(Box::new(BamSource {
        reader,
        header,
        text,
        record: bam::Record::default(),
        writer: sam::io::Writer::new(Vec::new()),
    }))
}
