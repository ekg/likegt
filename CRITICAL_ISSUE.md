# CRITICAL ISSUE: Graph Node ID Mismatch

## The Problem
The entire genotyping system is broken due to node ID mismatch between tools:

1. **odgi paths -H**: Outputs coverage for visited nodes only (1009 nodes)
   - Column headers are `node.1`, `node.2`, etc.
   - No information about which actual graph node IDs these correspond to
   - Nodes are likely in some internal order, not graph node ID order

2. **gafpack**: Outputs coverage for ALL nodes in the graph (3361 nodes)
   - Uses actual graph node IDs from the GFA
   - Includes nodes with 0 coverage

3. **Result**: Comparing coverage vectors that represent completely different nodes!
   - Node column 1 in odgi output ≠ Node 1 in gafpack output
   - This is why hold-2-out showed 25% identity (random)

## Why This Breaks Everything
- When we compare coverage vectors, we're comparing random nodes
- The cosine similarity is meaningless
- The genotyping results are essentially random

## Potential Solutions

### Option 1: Use consistent tools
- Use ONLY gafpack for both reference and test
- Generate reference by simulating reads from each haplotype
- Map all reads with the same pipeline

### Option 2: Use vg toolkit
- vg has consistent node handling
- Can extract coverage for all nodes consistently

### Option 3: Fix odgi output
- Need to determine which nodes odgi is outputting
- Map those to actual graph node IDs
- Align with gafpack output

## Current Status
**THE ENTIRE SYSTEM IS BROKEN**
- The Rust implementation assumes matching dimensions
- The Python scripts work around this incorrectly
- The genotyping results are meaningless

## Required Fix
Must ensure reference and test coverage vectors:
1. Have the same number of elements
2. Each element corresponds to the SAME graph node
3. Use consistent node ordering

Without this, the entire genotyping approach fails.