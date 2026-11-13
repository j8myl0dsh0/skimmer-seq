//! skimmer: streaming FASTQ quality filter.
//!
//! Pipeline per record: parse -> trim low-quality 3' tail -> length gate ->
//! mean-quality gate -> write. One record in memory at a time.

mod fastq;
mod quality;

use std::env;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process::ExitCode;

const USAGE: &str = "\
skimmer - streaming FASTQ quality filter

USAGE:
    skimmer [OPTIONS] [INPUT]

    INPUT defaults to stdin; records are written to stdout.

OPTIONS:
    --min-len <N>    drop records shorter than N bp after trimming  [default: 0]
    --min-mean <Q>   drop records with mean Phred below Q; 0 disables
                     -edge trimming                              [default: 0]
    --window <N>     window width for 3' tail trimming            [default: 20]
    -h, --help       print this message

EXIT CODES:
    0 success, 1 I/O or parse error, 2 bad usage";

/// How a run ended.
enum Failure {
    Usage(String),
    Runtime(String),
}

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Failure::Usage(m) | Failure::Runtime(m) => f.write_str(m),
        }
    }
}

#[derive(Default)]
struct Config {
    min_len: usize,
    min_mean: f64,
    window: usize,
}

#[derive(Default)]
struct Stats {
    read: u64,
    kept: u64,
    bases_in: u64,
    bases_out: u64,
    short: u64,
    low_quality: u64,
}

fn parse_args() -> Result<(Config, Option<String>), Failure> {
    let mut cfg = Config {
        min_len: 0,
        min_mean: 0.0,
        window: 20,
    };
    let mut input = None;
    let mut args = env::args().skip(1);

    while let Some(flag) = args.next() {
        let mut value = |name: &str| -> Result<String, Failure> {
            args.next()
                .ok_or_else(|| Failure::Usage(format!("{name} needs a value")))
        };

        match flag.as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                return Err(Failure::Usage(String::new()));
            }
            "--min-len" => {
                cfg.min_len = value("--min-len")?
                    .parse()
                    .map_err(|_| Failure::Usage("--min-len must be a whole number".into()))?;
            }
            "--min-mean" => {
                cfg.min_mean = value("--min-mean")?
                    .parse()
                    .map_err(|_| Failure::Usage("--min-mean must be a number".into()))?;
            }
            "--window" => {
                cfg.window = value("--window")?
                    .parse()
                    .map_err(|_| Failure::Usage("--window must be a whole number".into()))?;
                if cfg.window == 0 {
                    return Err(Failure::Usage("--window must be at least 1".into()));
                }
            }
            other if other.starts_with('-') && other != "-" => {
                return Err(Failure::Usage(format!("unknown option {other}")));
            }
            path => {
                if input.replace(path.to_string()).is_some() {
                    return Err(Failure::Usage("give at most one INPUT".into()));
                }
            }
        }
    }

    Ok((cfg, input))
}

fn run() -> Result<(), Failure> {
    let (cfg, input) = match parse_args() {
        Ok(parsed) => parsed,
        Err(Failure::Usage(m)) => {
            if m.is_empty() {
                return Ok(());
            }
            eprintln!("skimmer: {m}\n\n{USAGE}");
            return Err(Failure::Usage(m));
        }
        Err(other) => return Err(other),
    };

    let source: Box<dyn BufRead> = match input {
        Some(path) => {
            let file =
                File::open(&path).map_err(|e| Failure::Runtime(format!("open {path}: {e}")))?;
            Box::new(BufReader::with_capacity(1 << 20, file))
        }
        None => Box::new(BufReader::with_capacity(1 << 20, io::stdin())),
    };

    let mut records = fastq::Reader::new(source);
    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());

    let mut stats = Stats::default();
    while let Some(mut rec) = records
        .next_record()
        .map_err(|e| Failure::Runtime(e.to_string()))?
    {
        stats.read += 1;
        stats.bases_in += rec.seq.len() as u64;

        quality::trim_low_tail(&mut rec.seq, &mut rec.qual, cfg.window, cfg.min_mean);

        if rec.seq.len() < cfg.min_len {
            stats.short += 1;
            continue;
        }
        if quality::mean_phred(&rec.qual) < cfg.min_mean {
            stats.low_quality += 1;
            continue;
        }

        fastq::write_record(&mut out, &rec).map_err(|e| Failure::Runtime(e.to_string()))?;
        stats.kept += 1;
        stats.bases_out += rec.seq.len() as u64;
    }

    out.flush()
        .map_err(|e| Failure::Runtime(format!("write: {e}")))?;

    eprintln!(
        "read {} records ({} bp); kept {} ({} bp); dropped {} short, {} low-quality",
        stats.read, stats.bases_in, stats.kept, stats.bases_out, stats.short, stats.low_quality
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(Failure::Usage(_)) => ExitCode::from(2),
        Err(Failure::Runtime(m)) => {
            eprintln!("skimmer: {m}");
            ExitCode::from(1)
        }
    }
}
