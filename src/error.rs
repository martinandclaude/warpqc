//! Error type and the mapping to process exit codes (SPEC "Exit codes").

use std::fmt;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Bad header, coordinate-sorted or truncated input, reference mismatch,
    /// malformed record. Exit code 2.
    #[error("input error: {0}")]
    Input(String),
    /// I/O or resource failure. Exit code 3.
    #[error("I/O error: {context}: {source}")]
    Io {
        context: String,
        #[source]
        source: std::io::Error,
    },
}

impl Error {
    pub fn exit_code(&self) -> i32 {
        match self {
            Error::Input(_) => 2,
            Error::Io { .. } => 3,
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

pub fn input(msg: impl fmt::Display) -> Error {
    Error::Input(msg.to_string())
}

pub fn io(context: impl fmt::Display, source: std::io::Error) -> Error {
    Error::Io {
        context: context.to_string(),
        source,
    }
}

/// Attach a context string to an `io::Result`.
pub trait IoContext<T> {
    fn ctx(self, context: impl fmt::Display) -> Result<T>;
}

impl<T> IoContext<T> for std::io::Result<T> {
    fn ctx(self, context: impl fmt::Display) -> Result<T> {
        self.map_err(|e| io(context, e))
    }
}
