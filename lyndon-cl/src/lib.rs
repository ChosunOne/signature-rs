#[cfg(test)]
mod test {
    use std::marker::PhantomData;

    use cubecl::{
        benchmark::{Benchmark, TimingMethod},
        future,
        prelude::*,
        server::Handle,
        std::tensor::compact_strides,
        wgpu::{WgpuDevice, WgpuRuntime},
    };
    pub struct CpuTensor {
        pub data: Vec<f32>,
        pub strides: Vec<usize>,
        pub shape: Vec<usize>,
    }

    impl CpuTensor {
        pub fn arange(shape: Vec<usize>) -> Self {
            let size = shape.iter().product();
            let data = (0..size).map(|i| i as f32).collect();
            let strides = compact_strides(&shape);
            Self {
                data,
                strides,
                shape,
            }
        }

        pub fn empty(shape: Vec<usize>) -> Self {
            let size = shape.iter().product();
            let data = vec![0.0; size];
            let strides = compact_strides(&shape);
            Self {
                data,
                strides,
                shape,
            }
        }

        pub fn read(self) -> Vec<f32> {
            self.data
        }
    }

    pub fn reduce_matrix(input: &CpuTensor, output: &mut CpuTensor) {
        for i in 0..input.shape[0] {
            let mut acc = 0.0f32;
            for j in 0..input.shape[1] {
                acc += input.data[i * input.strides[0] + j];
            }
            output.data[i] = acc;
        }
    }

    #[derive(Debug)]
    pub struct GpuTensor<R: Runtime, F: Float + CubeElement> {
        data: Handle,
        shape: Vec<usize>,
        strides: Vec<usize>,
        _r: PhantomData<R>,
        _f: PhantomData<F>,
    }

    impl<R: Runtime, F: Float + CubeElement> Clone for GpuTensor<R, F> {
        fn clone(&self) -> Self {
            Self {
                data: self.data.clone(), // Handle is a pointer to the data, so cloning it is cheap
                shape: self.shape.clone(),
                strides: self.strides.clone(),
                _r: PhantomData,
                _f: PhantomData,
            }
        }
    }

    impl<R: Runtime, F: Float + CubeElement> GpuTensor<R, F> {
        /// Create a `GpuTensor` with a shape filled by number in order
        pub fn arange(shape: Vec<usize>, client: &ComputeClient<R::Server, R::Channel>) -> Self {
            let size = shape.iter().product();
            let data: Vec<F> = (0..size).map(|i| F::from_int(i as i64)).collect();
            let data = client.create(F::as_bytes(&data));

            let strides = cubecl::std::tensor::compact_strides(&shape);
            Self {
                data,
                shape,
                strides,
                _r: PhantomData,
                _f: PhantomData,
            }
        }

        /// Create an empty GpuTensor with a shape
        pub fn empty(shape: Vec<usize>, client: &ComputeClient<R::Server, R::Channel>) -> Self {
            let size = shape.iter().product::<usize>() * core::mem::size_of::<F>();
            let data = client.empty(size);

            let strides = compact_strides(&shape);
            Self {
                data,
                shape,
                strides,
                _r: PhantomData,
                _f: PhantomData,
            }
        }

        /// Create a TensorArg to pass to a kernel
        pub fn into_tensor_arg(&self, line_size: u8) -> TensorArg<'_, R> {
            unsafe {
                TensorArg::from_raw_parts::<F>(&self.data, &self.strides, &self.shape, line_size)
            }
        }

        /// Return the data from the client
        pub fn read(self, client: &ComputeClient<R::Server, R::Channel>) -> Vec<F> {
            let bytes = client.read_one(self.data.binding());
            F::from_bytes(&bytes).to_vec()
        }
    }

    pub struct ReductionBench<R: Runtime, F: Float + CubeElement> {
        input_shape: Vec<usize>,
        client: ComputeClient<R::Server, R::Channel>,
        _f: PhantomData<F>,
    }

    const LINE_SIZE: u32 = 4;

    impl<R: Runtime, F: Float + CubeElement> Benchmark for ReductionBench<R, F> {
        type Input = GpuTensor<R, F>;
        type Output = GpuTensor<R, F>;

        fn prepare(&self) -> Self::Input {
            GpuTensor::<R, F>::arange(self.input_shape.clone(), &self.client)
        }

        fn name(&self) -> String {
            format!("{}-reduction-{:?}", R::name(&self.client), self.input_shape).to_lowercase()
        }

        fn sync(&self) {
            future::block_on(self.client.sync());
        }

        fn execute(&self, input: Self::Input) -> Result<Self::Output, String> {
            let output_shape: Vec<usize> = vec![self.input_shape[0]];
            let output = GpuTensor::<R, F>::empty(output_shape, &self.client);

            unsafe {
                reduce_matrix_gpu::launch_unchecked::<F, R>(
                    &self.client,
                    CubeCount::Static(self.input_shape[0] as u32, 1, 1),
                    CubeDim::new(self.input_shape[1] as u32, 1, 1),
                    input.into_tensor_arg(LINE_SIZE as u8),
                    output.into_tensor_arg(LINE_SIZE as u8),
                );
            }

            Ok(output)
        }
    }

    #[cube(launch_unchecked)]
    pub fn reduce_matrix_gpu<F: Float>(input: &Tensor<Line<F>>, output: &mut Tensor<Line<F>>) {
        let mut acc = Line::new(F::new(0.0f32));
        for i in 0..input.shape(2) / LINE_SIZE {
            acc = acc + input[CUBE_POS_X * input.stride(0) + UNIT_POS_X * input.stride(1) + i];
        }
        output[CUBE_POS_X * output.stride(0) + UNIT_POS_X] = acc;
    }

    #[test]
    fn test_doing_stuff_on_the_cpu() {
        let input_shape = vec![3, 3];
        let output_shape = vec![3];
        let input = CpuTensor::arange(input_shape);
        let mut output = CpuTensor::empty(output_shape);

        reduce_matrix(&input, &mut output);

        println!("Executed reduction => {:?}", output.read());
    }

    #[test]
    fn test_doing_stuff_on_the_gpu() {
        let client = WgpuRuntime::client(&WgpuDevice::default());
        let input = GpuTensor::<WgpuRuntime, f32>::arange(vec![3, 3], &client);
        let output = GpuTensor::<WgpuRuntime, f32>::empty(vec![3, 3], &client);

        unsafe {
            reduce_matrix_gpu::launch_unchecked::<f32, WgpuRuntime>(
                &client,
                CubeCount::Static(1, 1, 1),
                CubeDim::new(1, 1, 1),
                input.into_tensor_arg(1),
                output.into_tensor_arg(1),
            );
        };

        println!(
            "Executed reduction with runtime {:?} => {:?}",
            WgpuRuntime::name(&client),
            output.read(&client)
        );
    }

    #[test]
    fn test_benchmark() {
        let client = WgpuRuntime::client(&WgpuDevice::default());
        let bench1 = ReductionBench::<WgpuRuntime, f32> {
            input_shape: vec![64, 256, 1024],
            client: client.clone(),
            _f: PhantomData,
        };

        let bench2 = ReductionBench::<WgpuRuntime, f32> {
            input_shape: vec![64, 64, 4096],
            client: client.clone(),
            _f: PhantomData,
        };

        for bench in [bench1, bench2] {
            println!("{}", bench.name());
            println!("{}", bench.run(TimingMethod::System).unwrap());
        }
    }
}
