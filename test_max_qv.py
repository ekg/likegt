#!/usr/bin/env python3
"""Quick test of max QV computation"""

import subprocess

def test_parse():
    # Get some PAF output
    cmd = ['allwave', '-i', 'hla-f.fa.gz', '-t', '4']
    result = subprocess.run(cmd, capture_output=True, text=True)
    paf_lines = result.stdout.strip().split('\n')
    
    print(f"Got {len(paf_lines)} PAF lines")
    
    # Look for HG00096 alignments
    hg00096_count = 0
    for line in paf_lines[:100]:  # Check first 100
        if line.startswith('HG00096#'):
            fields = line.split('\t')
            query = fields[0]
            target = fields[5]
            matches = int(fields[9])
            align_len = int(fields[10])
            identity = matches / align_len if align_len > 0 else 0
            
            query_id = query.split('#')[0] if '#' in query else query
            target_id = target.split('#')[0] if '#' in target else target
            
            if query_id == 'HG00096' and target_id != 'HG00096':
                hg00096_count += 1
                print(f"\nFound HG00096 alignment #{hg00096_count}:")
                print(f"  Query: {query[:50]}...")
                print(f"  Target: {target[:50]}...")
                print(f"  Identity: {identity:.4f}")
                print(f"  Matches: {matches}/{align_len}")
                
                # Compute QV
                error_rate = 1.0 - identity
                if error_rate <= 0:
                    qv = 60.0
                else:
                    import math
                    qv = -10.0 * math.log10(error_rate)
                    qv = min(60.0, max(0.0, qv))
                print(f"  QV: {qv:.1f}")
    
    print(f"\nTotal HG00096 non-self alignments found: {hg00096_count}")

if __name__ == '__main__':
    test_parse()