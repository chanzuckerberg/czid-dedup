use std::io::{Error, ErrorKind};

use super::fastx;

pub struct PairedRecords<T: fastx::Record, R: Iterator<Item = Result<T, std::io::Error>>> {
    records_r1: R,
    records_r2: R,
}

impl<T: fastx::Record, R: Iterator<Item = Result<T, std::io::Error>>> PairedRecords<T, R> {
    pub fn new(records_r1: R, records_r2: R) -> Self {
        PairedRecords {
            records_r1: records_r1,
            records_r2: records_r2,
        }
    }
}

impl<A: fastx::Record, T: Iterator<Item = Result<A, std::io::Error>>> Iterator
    for PairedRecords<A, T>
{
    type Item = Result<(A, A), Error>;

    fn next(&mut self) -> Option<Result<(A, A), Error>> {
        match (self.records_r1.next(), self.records_r2.next()) {
            (Some(Ok(r1_record)), Some(Ok(r2_record))) => {
                if r1_record.id() == r2_record.id() {
                    Some(Ok((r1_record, r2_record)))
                } else {
                    let message = format!(
                        "read pair had different read IDs: ({}, {})",
                        r1_record.id(),
                        r2_record.id()
                    );
                    Some(Err(Error::new(ErrorKind::InvalidData, message)))
                }
            }
            (None, None) => None,
            (Some(_), None) => Some(Err(Error::new(
                ErrorKind::UnexpectedEof,
                "reached the end of r2 before r1",
            ))),
            (None, Some(_)) => Some(Err(Error::new(
                ErrorKind::UnexpectedEof,
                "reached the end of r1 before r2",
            ))),
            (Some(Err(err)), _) => Some(Err(err)),
            (_, Some(Err(err))) => Some(Err(err)),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use bio::io::fasta;

    #[test]
    fn test_r1_longer() {
        let record = fasta::Record::with_attrs("id_a", None, &[]);
        let records_r1 = vec![Ok(record)].into_iter();
        let records_r2 = vec![].into_iter();
        let mut paired_iterator = PairedRecords::new(records_r1, records_r2);
        let result = paired_iterator.next();

        let error = result
            .expect("should return an element")
            .err()
            .expect("should return an error");
        assert_eq!(
            error.kind(),
            ErrorKind::UnexpectedEof,
            "should be of kind UnexpectedEof"
        );
        assert_eq!(
            error.to_string(),
            "reached the end of r2 before r1",
            "should contain correct message"
        );
    }

    #[test]
    fn test_r2_longer() {
        let record = fasta::Record::with_attrs("id_a", None, &[]);
        let records_r1 = vec![].into_iter();
        let records_r2 = vec![Ok(record)].into_iter();
        let mut paired_iterator = PairedRecords::new(records_r1, records_r2);
        let result = paired_iterator.next();

        let error = result
            .expect("should return an element")
            .err()
            .expect("should return an error");
        assert_eq!(
            error.kind(),
            ErrorKind::UnexpectedEof,
            "should be of kind UnexpectedEof"
        );
        assert_eq!(
            error.to_string(),
            "reached the end of r1 before r2",
            "should contain correct message"
        );
    }

    #[test]
    fn test_different_ids() {
        let record_r1 = fasta::Record::with_attrs("id_a", None, &[]);
        let record_r2 = fasta::Record::with_attrs("id_b", None, &[]);
        let records_r1 = vec![Ok(record_r1)].into_iter();
        let records_r2 = vec![Ok(record_r2)].into_iter();
        let mut paired_iterator = PairedRecords::new(records_r1, records_r2);
        let result = paired_iterator.next();

        let error = result
            .expect("should return an element")
            .err()
            .expect("should return an error");
        assert_eq!(
            error.kind(),
            ErrorKind::InvalidData,
            "should be of kind InvalidData"
        );
        assert_eq!(
            error.to_string(),
            "read pair had different read IDs: (id_a, id_b)",
            "should contain correct message"
        );
    }

    #[test]
    fn test_r1_error() {
        let records_r1 =
            vec![Err(Error::new(ErrorKind::Other, "I'm broken")) as Result<fasta::Record, Error>]
                .into_iter();
        let records_r2 = vec![Err(Error::new(ErrorKind::Other, "I'm also broken"))].into_iter();
        let mut paired_iterator = PairedRecords::new(records_r1, records_r2);
        let result = paired_iterator.next();

        let error = result
            .expect("should return an element")
            .err()
            .expect("should return an error");
        assert_eq!(error.kind(), ErrorKind::Other, "should be of kind Other");
        assert_eq!(error.to_string(), "I'm broken");
    }

    #[test]
    fn test_r2_error() {
        let record_r1 = fasta::Record::with_attrs("id_a", None, &[]);
        let records_r1 = vec![Ok(record_r1)].into_iter();
        let records_r2 = vec![Err(Error::new(ErrorKind::Other, "I'm broken"))].into_iter();
        let mut paired_iterator = PairedRecords::new(records_r1, records_r2);
        let result = paired_iterator.next();

        let error = result
            .expect("should return an element")
            .err()
            .expect("should return an error");
        assert_eq!(error.kind(), ErrorKind::Other, "should be of kind Other");
        assert_eq!(error.to_string(), "I'm broken");
    }
}
