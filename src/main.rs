use std::{
    error::Error,
    fmt::{Debug, Display},
    fs::File,
    ops::{AddAssign, MulAssign, SubAssign},
    path::{Path, PathBuf},
    str::FromStr,
};

use clap::Parser;
use file_type::FileType;
use ndarray::{Array2, s};
use num_traits::FromPrimitive;
use ordered_float::{FloatCore, NotNan, PrimitiveFloat};
use polars::{
    df,
    io::{SerReader, SerWriter},
    prelude::{
        CsvReadOptions, CsvWriter, Float32Type, Float64Type, IndexOrder, JsonFormat,
        JsonLineReader, JsonReader, JsonWriter, LazyFrame, ParquetWriter, PlPath,
        PolarsNumericType, ScanArgsParquet,
    },
};
use signature_rs::log_sig::LogSignatureBuilder;

/// Utility for loading path data from various file formats.
pub struct DataLoader;

impl DataLoader {
    /// Loads path data from a file, supporting multiple formats.
    ///
    /// Supports CSV, JSON, JSONL, and Parquet formats. The file format is
    /// automatically detected from the file's media type. Returns the path
    /// data as a 2D ndarray where rows are time points and columns are dimensions.
    pub fn load<T: PolarsNumericType>(path: &Path) -> Result<Array2<T::Native>, Box<dyn Error>> {
        let ft = *FileType::try_from_file(path)?
            .media_types()
            .first()
            .unwrap();
        let df = match ft {
            "text/csv" | "text/plain" => CsvReadOptions::default()
                .with_has_header(false)
                .try_into_reader_with_file_path(Some(path.into()))?
                .finish()?,
            "application/x-parquet" => LazyFrame::scan_parquet(
                PlPath::from_string(path.to_string_lossy().to_string()),
                ScanArgsParquet::default(),
            )?
            .collect()?,
            "application/json" => {
                let mut file = File::open(path)?;
                JsonReader::new(&mut file).finish()?
            }
            "application/x-ndjson"
            | "application/x-jsonlines"
            | "appliction/json-seq"
            | "application/jsonlines"
            | "application/jsonl"
            | "application/jsonlines+json" => JsonLineReader::from_path(path)?.finish()?,
            mt => return Err(format!("Unsupported file type: {mt}").into()),
        };

        Ok(df.to_ndarray::<T>(IndexOrder::Fortran)?)
    }
}

/// Supported output file formats for log signature results.
#[derive(Debug, Copy, Clone)]
pub enum OutputType {
    /// Standard JSON format.
    Json,
    /// JSON Lines format (one JSON object per line).
    JsonLines,
    /// Comma-separated values format.
    Csv,
    /// Apache Parquet columnar format.
    Parquet,
}

impl Display for OutputType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Json => write!(f, "json"),
            Self::JsonLines => write!(f, "jsonl"),
            Self::Csv => write!(f, "csv"),
            Self::Parquet => write!(f, "parquet"),
        }
    }
}

impl FromStr for OutputType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "json" => Ok(Self::Json),
            "jsonl" => Ok(Self::JsonLines),
            "csv" => Ok(Self::Csv),
            "parquet" => Ok(Self::Parquet),
            t => Err(format!("Unsupported output type: {t}")),
        }
    }
}

/// Supported numeric precision types for log signature computation.
#[derive(Debug, Copy, Clone)]
pub enum DataType {
    /// 32-bit floating-point precision.
    Float32,
    /// 64-bit floating-point precision.
    Float64,
}

impl FromStr for DataType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "f32" => Ok(Self::Float32),
            "f64" => Ok(Self::Float64),
            _ => Err(format!("Unsupported data type: {s}")),
        }
    }
}

impl Display for DataType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Float32 => write!(f, "f64"),
            Self::Float64 => write!(f, "f32"),
        }
    }
}

/// Command-line arguments for the log signature computation tool.
#[derive(Parser, Debug)]
pub struct Args {
    /// The number of dimensions of the path data.
    /// If not specified, will be inferred from the input data.
    #[arg(short = 'd')]
    num_dimensions: Option<usize>,
    /// The maximum degree for log signature computation.
    /// Higher degrees capture more complex geometric features.
    #[arg(short = 'k', default_value_t = 3)]
    max_degree: usize,
    /// The numeric precision to use for log signature coefficients.
    #[arg(short = 't', default_value_t = DataType::Float32)]
    data_type: DataType,
    /// The file path containing the input path data.
    #[arg(short = 'p')]
    path: PathBuf,
    /// The output file to write results to. If not specified, prints to stdout.
    #[arg(short = 'o')]
    output: Option<PathBuf>,
    /// The output file format.
    #[arg(short = 'f', default_value_t = OutputType::Csv)]
    output_type: OutputType,
}

/// Computes the log signature of a path using the specified parameters.
///
/// This function constructs a log signature builder with the given settings
/// and applies it to the input path data, returning the computed coefficients.
fn calculate_log_sig<
    T: AddAssign
        + Clone
        + Debug
        + Default
        + FloatCore
        + Send
        + Sync
        + FromPrimitive
        + PrimitiveFloat
        + SubAssign
        + MulAssign,
>(
    path: &Array2<T>,
    args: &Args,
) -> Result<Vec<NotNan<T>>, Box<dyn Error>> {
    let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));
    if path.shape().len() < 2 {
        return Err("Invalid path, expected at least two dimensions"
            .to_string()
            .into());
    }
    let num_dimensions = args.num_dimensions.unwrap_or(path.shape()[1]);
    let builder = LogSignatureBuilder::<u8>::new()
        .with_max_degree(args.max_degree)
        .with_num_dimensions(num_dimensions);
    let log_sig = builder.build_from_path(&path.slice(s![.., ..]));
    Ok(log_sig.series.coefficients)
}

#[allow(clippy::too_many_lines)]
pub fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    match args.data_type {
        DataType::Float32 => {
            let data = DataLoader::load::<Float32Type>(&args.path)?;
            let coefficients = calculate_log_sig(&data, &args)?
                .into_iter()
                .map(NotNan::into_inner)
                .collect::<Vec<_>>();

            if let Some(o) = args.output {
                let mut df = df!["coefficients" => coefficients]?;
                let mut file = File::create(&o)?;
                match args.output_type {
                    OutputType::Json => {
                        JsonWriter::new(&mut file).finish(&mut df)?;
                    }
                    OutputType::JsonLines => {
                        JsonWriter::new(&mut file)
                            .with_json_format(JsonFormat::JsonLines)
                            .finish(&mut df)?;
                    }
                    OutputType::Csv => {
                        CsvWriter::new(&mut file).finish(&mut df)?;
                    }
                    OutputType::Parquet => {
                        ParquetWriter::new(&mut file).finish(&mut df)?;
                    }
                }
            } else {
                println!("index \t coefficient");
                for (i, c) in coefficients.into_iter().enumerate() {
                    println!("{i}\t{c}");
                }
                println!();
            }
        }
        DataType::Float64 => {
            let data = DataLoader::load::<Float64Type>(&args.path)?;
            let coefficients = calculate_log_sig(&data, &args)?
                .into_iter()
                .map(NotNan::into_inner)
                .collect::<Vec<_>>();

            if let Some(o) = args.output {
                let mut df = df!["coefficients" => coefficients]?;
                let mut file = File::create(&o)?;
                match args.output_type {
                    OutputType::Json => {
                        JsonWriter::new(&mut file).finish(&mut df)?;
                    }
                    OutputType::JsonLines => {
                        JsonWriter::new(&mut file)
                            .with_json_format(JsonFormat::JsonLines)
                            .finish(&mut df)?;
                    }
                    OutputType::Csv => {
                        CsvWriter::new(&mut file).finish(&mut df)?;
                    }
                    OutputType::Parquet => {
                        ParquetWriter::new(&mut file).finish(&mut df)?;
                    }
                }
            } else {
                println!("index \t coefficient");
                for (i, c) in coefficients.into_iter().enumerate() {
                    println!("{i}\t{c}");
                }
                println!();
            }
        }
    }

    Ok(())
}
