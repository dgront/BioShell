use rustfft::{FftPlanner, num_complex::Complex};
use std::iter::zip;

pub fn autocorrelate_vectors(x: &Vec<f64>, y: &Vec<f64>, z: &Vec<f64>) -> Vec<f64> {
    let n = x.len();

    // Subtract the mean from the input data
    let mean_x = x.iter().sum::<f64>() / n as f64;
    let mean_y = y.iter().sum::<f64>() / n as f64;
    let mean_z = z.iter().sum::<f64>() / n as f64;

    let x: Vec<f64> = x.into_iter().map(|v| v - mean_x).collect();
    let y: Vec<f64> = y.into_iter().map(|v| v - mean_y).collect();
    let z: Vec<f64> = z.into_iter().map(|v| v - mean_z).collect();

    // Initialize planner for FFT
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(n);
    let ifft = planner.plan_fft_inverse(n);

    // Perform FFT on each component
    let mut x_fft: Vec<Complex<f64>> = x.into_iter().map(Complex::from).collect();
    let mut y_fft: Vec<Complex<f64>> = y.into_iter().map(Complex::from).collect();
    let mut z_fft: Vec<Complex<f64>> = z.into_iter().map(Complex::from).collect();

    fft.process(&mut x_fft);
    fft.process(&mut y_fft);
    fft.process(&mut z_fft);

    // Compute power spectrum (autocorrelation in frequency domain) for each axis
    let autocorr_x: Vec<Complex<f64>> = x_fft.iter().map(|&v| v * v.conj()).collect();
    let autocorr_y: Vec<Complex<f64>> = y_fft.iter().map(|&v| v * v.conj()).collect();
    let autocorr_z: Vec<Complex<f64>> = z_fft.iter().map(|&v| v * v.conj()).collect();

    // Perform inverse FFT to get the autocorrelation back in time domain
    let mut autocorr_x = autocorr_x;
    let mut autocorr_y = autocorr_y;
    let mut autocorr_z = autocorr_z;

    ifft.process(&mut autocorr_x);
    ifft.process(&mut autocorr_y);
    ifft.process(&mut autocorr_z);

    // Sum the autocorrelations of the three components
    let mut autocorr: Vec<f64> = zip(zip(autocorr_x.iter(), autocorr_y.iter()), autocorr_z.iter())
        .map(|((ax, ay), az)| (ax.re + ay.re + az.re) / n as f64)
        .collect();

    // Normalize and return only the real part of the autocorrelation
    let d = autocorr[0];
    autocorr.iter_mut().for_each(|v| *v /= d);

    return autocorr;
}

