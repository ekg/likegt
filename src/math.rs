
/// Calculate dot product of two vectors
pub fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Vectors must have same length");
    
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| x * y)
        .sum()
}

/// Calculate magnitude (Euclidean norm) of a vector
pub fn magnitude(v: &[f64]) -> f64 {
    v.iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
}

/// Calculate cosine similarity between two vectors
pub fn cosine_similarity(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Vectors must have same length");
    
    let dot = dot_product(a, b);
    let mag_a = magnitude(a);
    let mag_b = magnitude(b);
    
    if mag_a * mag_b > 0.0 {
        dot / (mag_a * mag_b)
    } else {
        0.0
    }
}

/// Sum multiple vectors element-wise
pub fn sum_vectors(vectors: &[&[f64]]) -> Vec<f64> {
    if vectors.is_empty() {
        return Vec::new();
    }
    
    let len = vectors[0].len();
    let mut result = vec![0.0; len];
    
    for vector in vectors {
        assert_eq!(vector.len(), len, "All vectors must have same length");
        for (i, &val) in vector.iter().enumerate() {
            result[i] += val;
        }
    }
    
    result
}

/// Calculate cosine similarity with optional masking
pub fn cosine_similarity_masked(a: &[f64], b: &[f64], mask: Option<&[bool]>) -> f64 {
    assert_eq!(a.len(), b.len(), "Vectors must have same length");
    
    if let Some(m) = mask {
        assert_eq!(a.len(), m.len(), "Mask must have same length as vectors");
        
        let mut dot = 0.0;
        let mut mag_a_sq = 0.0;
        let mut mag_b_sq = 0.0;
        
        for i in 0..a.len() {
            if m[i] {
                dot += a[i] * b[i];
                mag_a_sq += a[i] * a[i];
                mag_b_sq += b[i] * b[i];
            }
        }
        
        let mag = (mag_a_sq * mag_b_sq).sqrt();
        if mag > 0.0 {
            dot / mag
        } else {
            0.0
        }
    } else {
        cosine_similarity(a, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_dot_product() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        assert_eq!(dot_product(&a, &b), 32.0);
        
        // Test with zeros
        let c = vec![0.0, 0.0, 0.0];
        let d = vec![1.0, 2.0, 3.0];
        assert_eq!(dot_product(&c, &d), 0.0);
        
        // Test with negative values
        let e = vec![1.0, -2.0, 3.0];
        let f = vec![4.0, 5.0, -6.0];
        assert_eq!(dot_product(&e, &f), 1.0 * 4.0 + (-2.0) * 5.0 + 3.0 * (-6.0));
    }
    
    #[test]
    #[should_panic(expected = "Vectors must have same length")]
    fn test_dot_product_different_lengths() {
        let a = vec![1.0, 2.0];
        let b = vec![1.0, 2.0, 3.0];
        dot_product(&a, &b);
    }
    
    #[test]
    fn test_magnitude() {
        let v = vec![3.0, 4.0];
        assert_eq!(magnitude(&v), 5.0);
        
        // Test zero vector
        let zero = vec![0.0, 0.0, 0.0];
        assert_eq!(magnitude(&zero), 0.0);
        
        // Test unit vector
        let unit = vec![1.0, 0.0, 0.0];
        assert_eq!(magnitude(&unit), 1.0);
        
        // Test negative values
        let neg = vec![-3.0, -4.0];
        assert_eq!(magnitude(&neg), 5.0);
    }
    
    #[test]
    fn test_cosine_similarity() {
        // Orthogonal vectors
        let a = vec![1.0, 0.0];
        let b = vec![0.0, 1.0];
        assert_eq!(cosine_similarity(&a, &b), 0.0);
        
        // Identical vectors
        let c = vec![1.0, 1.0];
        let d = vec![1.0, 1.0];
        assert!((cosine_similarity(&c, &d) - 1.0).abs() < 1e-10);
        
        // Opposite vectors
        let e = vec![1.0, 0.0];
        let f = vec![-1.0, 0.0];
        assert!((cosine_similarity(&e, &f) + 1.0).abs() < 1e-10);
        
        // Zero vectors
        let g = vec![0.0, 0.0];
        let h = vec![1.0, 1.0];
        assert_eq!(cosine_similarity(&g, &h), 0.0);
        
        // Real-world example
        let v1 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let v2 = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((cosine_similarity(&v1, &v2) - 1.0).abs() < 1e-10); // Parallel vectors
    }
    
    #[test]
    fn test_sum_vectors() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let result = sum_vectors(&[&a, &b]);
        assert_eq!(result, vec![5.0, 7.0, 9.0]);
        
        // Test with single vector
        let single = sum_vectors(&[&a]);
        assert_eq!(single, vec![1.0, 2.0, 3.0]);
        
        // Test with three vectors
        let c = vec![7.0, 8.0, 9.0];
        let triple = sum_vectors(&[&a, &b, &c]);
        assert_eq!(triple, vec![12.0, 15.0, 18.0]);
        
        // Test empty input
        let empty = sum_vectors(&[]);
        assert_eq!(empty, Vec::<f64>::new());
    }
    
    #[test]
    fn test_cosine_similarity_masked() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![2.0, 4.0, 6.0, 8.0];
        
        // No mask - should be identical
        assert_eq!(
            cosine_similarity_masked(&a, &b, None),
            cosine_similarity(&a, &b)
        );
        
        // Mask out middle elements
        let mask = vec![true, false, false, true];
        let similarity = cosine_similarity_masked(&a, &b, Some(&mask));
        // Should only consider elements 0 and 3
        let expected = (1.0 * 2.0 + 4.0 * 8.0) / ((1.0_f64 * 1.0 + 4.0 * 4.0).sqrt() * (2.0_f64 * 2.0 + 8.0 * 8.0).sqrt());
        assert!((similarity - expected).abs() < 1e-10);
        
        // All masked out
        let all_false = vec![false, false, false, false];
        assert_eq!(cosine_similarity_masked(&a, &b, Some(&all_false)), 0.0);
    }
    
    #[test]
    #[should_panic(expected = "Mask must have same length")]
    fn test_cosine_similarity_masked_wrong_length() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let mask = vec![true, false];
        cosine_similarity_masked(&a, &b, Some(&mask));
    }
}