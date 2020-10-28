use bio::io::{fasta, fastq};
use clap::{App, Arg};
use simple_error;
use std::error::Error;
use std::fs::File;

mod clusters;
mod fastx;
mod paired;

macro_rules! box_result_error {
    ($result:expr) => {
        $result.map_err(Box::new)
    };
}

macro_rules! unwrap_or_return {
    ($result:expr) => {{
        match $result {
            Err(err) => return Err(err),
            Ok(v) => v,
        }
    }};
}

macro_rules! box_bail {
    ($result:expr) => {
        unwrap_or_return!(box_result_error!($result))
    };
}

macro_rules! dedup {
    ($fastx:tt, $fastx_type_r1:expr, $input_r1:expr, $output_r1:expr, $inputs:expr, $outputs:expr, $clusters:expr) => {{
        let records_r1 = $fastx::Reader::from_file($input_r1).unwrap().records();
        let writer_r1 = $fastx::Writer::to_file($output_r1).unwrap();
        match ($inputs.next(), $outputs.next()) {
            (Some(input_r2), Some(output_r2)) => {
                let fastx_type_r2 = fastx::fastx_type(input_r2).unwrap();
                if fastx_type_r2 != $fastx_type_r1 {
                    let message = format!(
                        "paired inputs have different file types r1: {}, r2: {}",
                        $fastx_type_r1, fastx_type_r2
                    );
                    return Err(Box::new(simple_error::simple_error!(message)));
                }
                let records_r2 = $fastx::Reader::from_file(input_r2).unwrap().records();
                let writer_r2 = $fastx::Writer::to_file(output_r2).unwrap();
                let records = paired::PairedRecords::new(records_r1, records_r2);
                pair(records, writer_r1, writer_r2, $clusters)
            }
            (None, None) => single(records_r1, writer_r1, $clusters),
            _ => panic!("must have the same number of inputs and outputs"),
        }
    }};
}

fn single<
    T: fastx::Record,
    R: Iterator<Item = Result<T, std::io::Error>>,
    S: fastx::Writer<T>,
    U: std::io::Write,
>(
    records: R,
    mut writer: S,
    mut clusters: clusters::Clusters<U>,
) -> Result<clusters::Clusters<U>, Box<dyn Error>> {
    for result in records {
        let record = box_bail!(result);
        box_bail!(record
            .check()
            .map_err(|err| simple_error::simple_error!(err)));

        let result = clusters.insert_single(&record);
        if box_bail!(result) {
            box_bail!(writer.write_record(&record));
        }
    }
    Ok(clusters)
}

fn pair<
    T: fastx::Record,
    R: Iterator<Item = Result<T, std::io::Error>>,
    S: fastx::Writer<T>,
    U: std::io::Write,
>(
    records: paired::PairedRecords<T, R>,
    mut writer_r1: S,
    mut writer_r2: S,
    mut clusters: clusters::Clusters<U>,
) -> Result<clusters::Clusters<U>, Box<dyn Error>> {
    for result in records {
        let record = box_bail!(result);

        box_bail!(record
            .check()
            .map_err(|err| simple_error::simple_error!(&err)));

        let result = clusters.insert_pair(&record);
        if box_bail!(result) {
            box_bail!(writer_r1.write_record(record.r1()));
            box_bail!(writer_r2.write_record(record.r2()));
        }
    }
    Ok(clusters)
}

fn run_dedup<T: Into<std::ffi::OsString> + Clone, R: IntoIterator<Item = T>>(
    args: R,
) -> Result<clusters::Clusters<File>, Box<dyn Error>> {
    let matches = App::new(clap::crate_name!())
        .version(clap::crate_version!())
        .author(clap::crate_authors!())
        .about(clap::crate_description!())
        .arg(
            Arg::with_name("inputs")
                .short("i")
                .long("inputs")
                .help("Input FASTQ")
                .multiple(true)
                .min_values(1)
                .max_values(2)
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("deduped-outputs")
                .short("o")
                .long("deduped-outputs")
                .help("Output deduped FASTQ")
                .multiple(true)
                .min_values(1)
                .max_values(2)
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("cluster-output")
                .short("c")
                .long("cluster-output")
                .help("Output cluster file")
                .takes_value(true)
                .default_value("clusters.csv"),
        )
        .arg(
            Arg::with_name("prefix-length")
                .short("l")
                .long("prefix-length")
                .help("Length of the prefix to consider")
                .takes_value(true),
        )
        .get_matches_from(args);

    // presence guarunteed by clap
    let mut inputs = matches.values_of("inputs").unwrap();
    let mut outputs = matches.values_of("deduped-outputs").unwrap();
    let cluster_output = matches.value_of("cluster-output").unwrap();
    let prefix_length_opt = matches
        .value_of("prefix-length")
        .map(|n| n.parse::<usize>().unwrap());
    let input_r1 = inputs.next().unwrap();
    let output_r1 = outputs.next().unwrap();

    let bytes = File::open(input_r1).unwrap().metadata().unwrap().len() as usize;
    // 400 is based on the bytes per record of an example file, should be reasonable
    let clusters =
        clusters::Clusters::from_file(cluster_output, prefix_length_opt, bytes / 400).unwrap();

    match fastx::fastx_type(input_r1).unwrap() {
        fastx::FastxType::Fasta => dedup!(
            fasta,
            fastx::FastxType::Fasta,
            input_r1,
            output_r1,
            inputs,
            outputs,
            clusters
        ),
        fastx::FastxType::Fastq => dedup!(
            fastq,
            fastx::FastxType::Fastq,
            input_r1,
            output_r1,
            inputs,
            outputs,
            clusters
        ),
        fastx::FastxType::Invalid => Err(Box::new(simple_error::simple_error!(
            "input file is not a valid FASTA or FASTQ file"
        )) as Box<dyn Error>),
    }
}

fn main() {
    match run_dedup(std::env::args()) {
        Err(err) => println!("{}", err.to_string()),
        Ok(info) => {
            println!(
                "duplicates:   {:width$}",
                info.duplicate_records(),
                width = 16
            );
            println!("unique reads: {:width$}", info.unique_records(), width = 16);
            println!("total reads:  {:width$}", info.total_records(), width = 16);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use bio::io::fastq;
    use rand::Rng;
    use tempfile::tempdir;

    fn random_seq(len: usize) -> Vec<u8> {
        const CHARSET: &[u8] = b"ACTG";
        let mut rng = rand::thread_rng();
        (0..len)
            .map(|_| {
                let idx = rng.gen_range(0, CHARSET.len());
                CHARSET[idx]
            })
            .collect()
    }

    #[test]
    fn test_run_dedup_single() {
        let dir = tempdir().unwrap();
        let input_path = dir.path().join("input.fastq").to_str().unwrap().to_string();
        let output_path = dir
            .path()
            .join("output.fastq")
            .to_str()
            .unwrap()
            .to_string();
        let cluster_path = dir.path().join("cluster.csv").to_str().unwrap().to_string();

        {
            let mut writer = fastq::Writer::to_file(&input_path).expect("don't break");
            let seq = random_seq(20);
            writer.write("id_a", None, &seq, &seq).expect("don't break");
        }

        let args = [
            "executable",
            "-i",
            &input_path,
            "-o",
            &output_path,
            "-c",
            &cluster_path,
        ];
        let result = run_dedup(&args).expect("don't break");
        assert_eq!(result.total_records(), 1);
        dir.close().expect("don't break");
    }

    #[test]
    fn test_run_dedup_paired() {
        let dir = tempdir().unwrap();
        let input_path_r1 = dir
            .path()
            .join("input-r1.fasta")
            .to_str()
            .unwrap()
            .to_string();
        let input_path_r2 = dir
            .path()
            .join("input-r2.fasta")
            .to_str()
            .unwrap()
            .to_string();
        let output_path_r1 = dir
            .path()
            .join("output-r1.fasta")
            .to_str()
            .unwrap()
            .to_string();
        let output_path_r2 = dir
            .path()
            .join("output-r2.fasta")
            .to_str()
            .unwrap()
            .to_string();
        let cluster_path = dir.path().join("cluster.csv").to_str().unwrap().to_string();

        {
            let mut writer_r1 = fasta::Writer::to_file(&input_path_r1).expect("don't break");
            let mut writer_r2 = fasta::Writer::to_file(&input_path_r2).expect("don't break");
            let seq = random_seq(20);
            writer_r1.write("id_a", None, &seq).expect("don't break");
            writer_r2.write("id_a", None, &seq).expect("don't break");
        }

        let args = [
            "executable",
            "-i",
            &input_path_r1,
            "-i",
            &input_path_r2,
            "-o",
            &output_path_r1,
            "-o",
            &output_path_r2,
            "-c",
            &cluster_path,
        ];
        let result = run_dedup(&args).expect("don't break");
        assert_eq!(result.total_records(), 1);
        dir.close().expect("don't break");
    }

    #[test]
    fn test_run_dedup_paired_mismatched_files() {
        let dir = tempdir().unwrap();
        let input_path_r1 = dir
            .path()
            .join("input-r1.fasta")
            .to_str()
            .unwrap()
            .to_string();
        let input_path_r2 = dir
            .path()
            .join("input-r2.fastq")
            .to_str()
            .unwrap()
            .to_string();
        let output_path_r1 = dir
            .path()
            .join("output-r1.fasta")
            .to_str()
            .unwrap()
            .to_string();
        let output_path_r2 = dir
            .path()
            .join("output-r2.fasta")
            .to_str()
            .unwrap()
            .to_string();
        let cluster_path = dir.path().join("cluster.csv").to_str().unwrap().to_string();

        {
            let mut writer_r1 = fasta::Writer::to_file(&input_path_r1).expect("don't break");
            let mut writer_r2 = fastq::Writer::to_file(&input_path_r2).expect("don't break");
            let seq = random_seq(20);
            writer_r1.write("id_a", None, &seq).expect("don't break");
            writer_r2
                .write("id_a", None, &seq, &seq)
                .expect("don't break");
        }

        let args = [
            "executable",
            "-i",
            &input_path_r1,
            "-i",
            &input_path_r2,
            "-o",
            &output_path_r1,
            "-o",
            &output_path_r2,
            "-c",
            &cluster_path,
        ];
        let result = run_dedup(&args);
        let message = result
            .err()
            .expect("should error on mismatched inputs")
            .to_string();
        assert_eq!(
            message,
            "paired inputs have different file types r1: fasta, r2: fastq"
        );
        dir.close().expect("don't break");
    }
}
