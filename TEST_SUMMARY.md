# Test Suite Summary for LikeGT

## Overview
Comprehensive test suite with **20 tests** covering all major components of the genotyping system.

## Test Categories

### 1. Unit Tests

#### Math Module (`src/math.rs`)
- ✅ `test_dot_product`: Validates dot product calculation
- ✅ `test_dot_product_different_lengths`: Ensures panic on mismatched vector lengths
- ✅ `test_magnitude`: Tests Euclidean norm calculation
- ✅ `test_cosine_similarity`: Validates cosine similarity for various vector pairs
- ✅ `test_sum_vectors`: Tests vector addition
- ✅ `test_cosine_similarity_masked`: Tests masked cosine similarity
- ✅ `test_cosine_similarity_masked_wrong_length`: Ensures panic on wrong mask length

#### IO Module (`src/io.rs`)
- ✅ `test_coverage_data_creation`: Tests CoverageData struct creation

#### Coverage Module (`src/coverage.rs`)
- ✅ `test_genotype_data_creation`: Tests GenotypeData initialization
- ✅ `test_genotype_data_mismatched_lengths`: Validates error on mismatched coverage lengths
- ✅ `test_filter_blacklist`: Tests filtering haplotypes by blacklist
- ✅ `test_empty_blacklist`: Tests behavior with empty blacklist

#### Genotype Module (`src/genotype.rs`)
- ✅ `test_combinations_count`: Validates combination counting formula

### 2. Integration Tests (`tests/integration_test.rs`)
- ✅ `test_integration_small_example`: End-to-end test with 5 haplotypes
- ✅ `test_perfect_genotyping_scenario`: Tests multiple genotyping scenarios
- ✅ `test_read_test_data_if_exists`: Tests reading real HLA data

### 3. Hold-0-Out Validation (`tests/hold0_validation_test.rs`)
- ✅ `test_hold0_always_correct`: Validates perfect recovery with exact coverage
- ✅ `test_hold0_with_small_noise`: Tests robustness to small noise

### 4. Library Tests (`src/lib.rs`)
- ✅ `test_hold0_perfect_recovery`: Tests perfect genotype recovery
- ✅ `test_coverage_scaling_invariance`: Validates cosine similarity scale invariance

## Key Test Properties

### Correctness Guarantees
1. **Perfect Recovery**: When given exact coverage data (hold-0-out), the system always recovers the correct genotype with similarity = 1.0
2. **Scale Invariance**: Cosine similarity is invariant to coverage scaling
3. **Noise Robustness**: Small noise (< 5%) still allows correct genotype recovery

### Edge Cases Covered
- Empty inputs
- Mismatched vector lengths
- Homozygous genotypes
- Zero vectors
- Masked similarity calculations

## Running Tests

```bash
# Run all tests
cargo test

# Run specific test module
cargo test math::tests

# Run integration tests only
cargo test --test integration_test

# Run with output
cargo test -- --nocapture

# Run in release mode (faster)
cargo test --release
```

## Test Data
- Synthetic data generated programmatically for unit tests
- Real HLA-F graph data in `tests/data/` for integration testing
- Test coverage matrices with known ground truth

## Coverage Validation
The test suite validates that:
1. Mathematical operations are correct
2. File I/O handles gzip TSV correctly
3. Genotyping algorithm recovers true genotypes
4. System is robust to realistic data variations
5. Hold-0-out always produces perfect results with exact data