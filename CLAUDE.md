# Claude Development Guide

## Journaling

When debugging complex issues or implementing algorithms, I maintain a development journal at `JOURNAL.md`. This helps me:

1. **Think through problems** without cluttering the conversation
2. **Track what's actually broken** vs what appears broken  
3. **Remember the real goal** instead of getting lost in workarounds
4. **Reflect on mistakes** and course-correct

The journal is for my own reasoning process. You don't need to read it unless curious.

## Running Tests

### Hold-2-out validation
Tests genotyping accuracy by holding out one individual at a time:

```bash
# Single individual
cargo run -- hold2-out -f hla-f.fa.gz -g hla-f.k51.gfa -i HG00096

# With verbose output
cargo run -- hold2-out -f hla-f.fa.gz -g hla-f.k51.gfa -i HG00096 --verbose
```

## Key Commands

- `cargo test` - Run unit tests
- `cargo run -- --help` - See all subcommands
- `cargo build --release` - Build optimized binary

## Architecture Notes

The toolkit implements the COSIGT algorithm for pangenome genotyping. Core components:

- **io.rs**: Reading coverage matrices
- **validation.rs**: Hold-k-out validation logic  
- **math.rs**: Cosine similarity calculations
- **commands/hold2out.rs**: Full pipeline for validation

## Current Focus

Implementing proper COSIGT validation with:
- Real node coverage from aligned reads
- Correct graph coordinate mapping
- Batch processing capabilities