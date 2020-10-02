use core::hash::Hash;
use core::hash::Hasher;
use std::cmp;
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::io;

use super::fastx;

pub struct Clusters<T: io::Write> {
    cluster_map: HashMap<u64, String>, 
    cluster_csv_writer: csv::Writer<T>,
    total_records: u64,
    prefix_length_opt: Option<usize>,
}

impl<T: std::io::Write> Clusters<T> {
    fn insert_record(&mut self, seq_hash: u64, id: String) -> Result<bool, csv::Error> {
        self.total_records += 1;
        match self.cluster_map.get(&seq_hash) {
            Some(existing_id) => self.cluster_csv_writer.write_record(vec![existing_id, &id]).map(|_| false),
            None => {
                let res = self.cluster_csv_writer.write_record(vec![&id, &id]);
                self.cluster_map.insert(seq_hash, id);
                res.map(|_| true)
            }
        }
    }

    fn get_prefix<'a, 'b>(&'a self, seq: &'b[u8]) -> &'b[u8] {
        let seq_length = seq.len();
        let prefix_length = self.prefix_length_opt
            .map(|prefix_length| cmp::min(prefix_length, seq_length))
            .unwrap_or(seq_length);
        &seq[..prefix_length]
    }

    pub fn insert_single<R: fastx::Record>(&mut self, record: &R) -> Result<bool, csv::Error> {
        let mut seq_hasher = DefaultHasher::new();
        Hash::hash_slice(self.get_prefix(record.seq()), &mut seq_hasher);
        let seq_hash = seq_hasher.finish();
        self.insert_record(seq_hash, record.id().to_owned())
    }

    pub fn insert_pair<R: fastx::Record>(&mut self, record_r1: &R, record_r2: &R) -> Result<bool, csv::Error> {
        let mut seq_hasher = DefaultHasher::new();
        Hash::hash_slice(self.get_prefix(record_r1.seq()), &mut seq_hasher);
        Hash::hash(&0, &mut seq_hasher);
        Hash::hash_slice(self.get_prefix(record_r2.seq()), &mut seq_hasher);
        let seq_hash = seq_hasher.finish();
        self.insert_record(seq_hash, record_r1.id().to_owned())
    }

    pub fn unique_records(&self) -> u64 {
        self.cluster_map.len() as u64
    }

    pub fn duplicate_records(&self) -> u64 {
        self.total_records - self.unique_records()
    }

    pub fn total_records(&self) -> u64 {
        self.total_records
    }

    pub fn from_writer(cluster_output: T, prefix_length_opt: Option<usize>, capacity: usize) -> Result<Self, csv::Error> {
        let mut cluster_csv_writer = csv::Writer::from_writer(cluster_output);
        cluster_csv_writer.write_record(vec!["representative read id", "read id"]).map(|_| {
            let cluster_map = HashMap::with_capacity(capacity);
            Clusters { cluster_map: cluster_map, cluster_csv_writer: cluster_csv_writer, total_records: 0, prefix_length_opt: prefix_length_opt }
        })
    }
}

impl Clusters<File> {
    pub fn from_file<P: AsRef<std::path::Path>>(cluster_output_path: P, prefix_length_opt: Option<usize>, capacity: usize) -> Result<Self, csv::Error> {
        File::create(cluster_output_path)
            .map_err(csv::Error::from)
            .and_then(|cluster_output| Clusters::from_writer(cluster_output, prefix_length_opt, capacity))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use bio::io::fasta;
    use std::io::Cursor;
    use std::str;
    use rand::Rng;

    fn random_seq(len: usize) -> Vec<u8> {
        const CHARSET: &[u8] = b"ACTG";
        let mut rng = rand::thread_rng();
        (0..len)
        .map(|_| {
            let idx = rng.gen_range(0, CHARSET.len());
            CHARSET[idx]
        }).collect()
    }

    #[test]
    fn test_insert_single() {
        let mut cluster_output = Cursor::new(Vec::new());
        {
            let mut clusters = Clusters::from_writer(&mut cluster_output, Some(10), 200).expect("asdasd");
            let seq = random_seq(20);
            let record_1 = fasta::Record::with_attrs("id_a", None, &seq);
            clusters.insert_single(&record_1).expect("don't break");
            let record_2 = fasta::Record::with_attrs("id_b", None, &seq);
            clusters.insert_single(&record_2).expect("don't break");
            assert_eq!(clusters.duplicate_records(), 1);
            assert_eq!(clusters.unique_records(), 1);
            assert_eq!(clusters.total_records(), 2);
        }
        assert_eq!(str::from_utf8(cluster_output.into_inner().as_slice()).unwrap(), "representative read id,read id\nid_a,id_a\nid_a,id_b\n");
    }

    #[test]
    fn test_insert_pair() {
        let mut cluster_output = Cursor::new(Vec::new());
        {
            let mut clusters = Clusters::from_writer(&mut cluster_output, Some(10), 200).expect("asdasd");
            let seq_r1 = random_seq(20);
            let seq_r2 = random_seq(20);
            let record_1_r1 = fasta::Record::with_attrs("id_a", None, &seq_r1);
            let record_1_r2 = fasta::Record::with_attrs("id_a", None, &seq_r2);
            clusters.insert_pair(&record_1_r1, &record_1_r2).expect("don't break");
            let record_2_r1 = fasta::Record::with_attrs("id_b", None, &seq_r1);
            let record_2_r2 = fasta::Record::with_attrs("id_b", None, &seq_r2);
            clusters.insert_pair(&record_2_r1, &record_2_r2).expect("don't break");
            assert_eq!(clusters.duplicate_records(), 1);
            assert_eq!(clusters.unique_records(), 1);
            assert_eq!(clusters.total_records(), 2);
        }
        assert_eq!(str::from_utf8(cluster_output.into_inner().as_slice()).unwrap(), "representative read id,read id\nid_a,id_a\nid_a,id_b\n");
    }


}
