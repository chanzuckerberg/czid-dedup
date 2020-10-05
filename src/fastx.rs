use bio::io::{fasta, fastq};
use std::fs::File;
use std::io::prelude::*;
use std::io::Write;

pub trait Record {
    fn id(&self) -> &str;
    fn seq(&self) -> &[u8];
    fn check(&self) -> Result<(), &str>;
}

impl Record for fasta::Record {
    fn id(&self) -> &str {
        self.id()
    }

    fn seq(&self) -> &[u8] {
        self.seq()
    }

    fn check(&self) -> Result<(), &str> {
        self.check()
    }
}

impl Record for fastq::Record {
    fn id(&self) -> &str {
        self.id()
    }

    fn seq(&self) -> &[u8] {
        self.seq()
    }

    fn check(&self) -> Result<(), &str> {
        self.check()
    }
}

pub trait Writer<T: Record> {
    fn write_record(&mut self, record: &T) -> Result<(), std::io::Error>;
}

impl<T: Write> Writer<fasta::Record> for fasta::Writer<T> {
    fn write_record(&mut self, record: &fasta::Record) -> Result<(), std::io::Error> {
        self.write_record(&record)
    }
}

impl<T: Write> Writer<fastq::Record> for fastq::Writer<T> {
    fn write_record(&mut self, record: &fastq::Record) -> Result<(), std::io::Error> {
        self.write_record(&record)
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum FastxType {
    Fastq,
    Fasta,
    Invalid,
}

pub fn fastx_type<P: AsRef<std::path::Path>>(path: P) -> Result<FastxType, std::io::Error> {
    let mut file = match File::open(path) {
        Ok(f) => f,
        Err(err) => return Err(err),
    };
    let mut byte = [0; 1];
    if let Err(err) = file.read(&mut byte) {
        return Err(err);
    }

    match byte[0] as char {
        '>' => Ok(FastxType::Fasta),
        '@' => Ok(FastxType::Fastq),
        _ => Ok(FastxType::Invalid),
    }
}

impl std::fmt::Display for FastxType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            FastxType::Fasta => "fasta",
            FastxType::Fastq => "fastq",
            FastxType::Invalid => "invalid",
        };
        write!(f, "{}", s)
    }
}
