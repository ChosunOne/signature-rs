use cubecl::{prelude::*, std::tensor::TensorHandle};
use ordered_float::NotNan;
use std::marker::PhantomData;

use signature_rs::{LogSignature, LogSignatureBuilder};

use crate::log_signature::LogSignatureCLData;

pub struct LogSignatureCLBuilder<R> {
    log_signature_structure: LogSignature<u8, NotNan<f32>>,
    _r: PhantomData<R>,
}

impl<R: Runtime> LogSignatureCLBuilder<R> {
    #[must_use]
    pub fn new(max_degree: usize, num_dimensions: usize) -> Self {
        let builder = LogSignatureBuilder::<u8>::new()
            .with_max_degree(max_degree)
            .with_num_dimensions(num_dimensions);
        let log_signature_structure = builder.build::<NotNan<f32>>();
        Self {
            log_signature_structure,
            _r: PhantomData,
        }
    }

    pub fn build_from_tensor_path<'a>(
        &'a self,
        path: &TensorHandle<R, f32>,
    ) -> LogSignatureCLData<'a, R, f32, NotNan<f32>> {
        todo!()
    }
}
