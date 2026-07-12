use pyo3::prelude::*;

/// Euclidean (L2) norm of a vector.
#[pyfunction]
fn l2_norm(v: Vec<f64>) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[pymodule]
fn qscat_kernels(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(l2_norm, m)?)?;
    Ok(())
}
