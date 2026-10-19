use pyo3::exceptions::{PyMemoryError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

unsafe extern "C" {
    fn qscat_ld_size() -> usize;
    fn qscat_ld_bits() -> u32;
    fn qscat_ld_exponent() -> u32;
    fn qscat_ld_supported() -> i32;
    fn qscat_ld_residual(
        a: *const u8,
        b: *const u8,
        x: *const u8,
        out: *mut u8,
        m: usize,
        k: usize,
        n: usize,
    ) -> i32;
    fn qscat_ld_fixture(out: *mut u8, entries: usize);
}

/// Native scalar storage, significand, exponent range and supported ABI.
#[pyfunction]
pub fn capability() -> (usize, u32, u32, bool) {
    // These C functions return compile-time scalar properties and access no pointers.
    unsafe {
        (
            qscat_ld_size(),
            qscat_ld_bits(),
            qscat_ld_exponent(),
            qscat_ld_supported() != 0,
        )
    }
}

fn byte_count(rows: usize, cols: usize) -> Result<usize, &'static str> {
    if rows == 0 || cols == 0 {
        return Err("Positive dimensions required");
    }
    rows.checked_mul(cols)
        .and_then(|n| n.checked_mul(capability().0))
        .filter(|&n| n <= isize::MAX as usize)
        .ok_or("Matrix byte count overflow")
}

/// Internal checked byte ABI; Python numerical interfaces retain array signatures.
pub fn checked_residual(
    a: &[u8],
    b: &[u8],
    x: &[u8],
    m: usize,
    k: usize,
    n: usize,
) -> Result<Vec<u8>, &'static str> {
    if !capability().3 {
        return Err("Unsupported native long-double format");
    }
    let a_size = byte_count(m, k)?;
    let b_size = byte_count(m, n)?;
    let x_size = byte_count(k, n)?;
    if (a.len(), b.len(), x.len()) != (a_size, b_size, x_size) {
        return Err("Matrix buffer length differs from dimensions");
    }
    let mut output = Vec::new();
    output
        .try_reserve_exact(b_size)
        .map_err(|_| "Allocation failed")?;
    output.resize(b_size, 0);
    // Checked sizes prove every C input/output range and allocation product valid.
    // C uses aligned internal copies and memcpy for the possibly unaligned buffers.
    let status = unsafe {
        qscat_ld_residual(
            a.as_ptr(),
            b.as_ptr(),
            x.as_ptr(),
            output.as_mut_ptr(),
            m,
            k,
            n,
        )
    };
    if status != 0 {
        return Err("Native residual allocation or format failure");
    }
    Ok(output)
}

#[pyfunction]
fn residual_bytes<'py>(
    py: Python<'py>,
    a: &[u8],
    b: &[u8],
    x: &[u8],
    m: usize,
    k: usize,
    n: usize,
) -> PyResult<Bound<'py, PyBytes>> {
    let result = py
        .detach(|| checked_residual(a, b, x, m, k, n))
        .map_err(|error| {
            if error.contains("Allocation") {
                PyMemoryError::new_err(error)
            } else {
                PyValueError::new_err(error)
            }
        })?;
    Ok(PyBytes::new(py, &result))
}

/// Exact dyadic dense inputs for Criterion; not a scientific validation oracle.
pub fn benchmark_fixture(entries: usize) -> Vec<u8> {
    assert!(capability().3);
    let mut result = vec![0; byte_count(entries, 1).expect("Fixture size")];
    // The destination has the checked space for exactly entries scalars.
    unsafe { qscat_ld_fixture(result.as_mut_ptr(), entries) };
    result
}

#[pymodule]
fn qscat_extended_kernels(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(capability, m)?)?;
    m.add_function(wrap_pyfunction!(residual_bytes, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn checked_dimensions_reject_empty_and_overflow() {
        assert!(byte_count(0, 1).is_err());
        assert!(byte_count(usize::MAX, 2).is_err());
    }
    #[test]
    fn short_buffers_are_rejected() {
        assert!(checked_residual(&[], &[], &[], 1, 1, 1).is_err());
    }
}
