import nupack
import os
import sys

# NUPACK inspection
test_structure = '((((....))))'
config = nupack.Model(material='rna')
domain = nupack.Domain('N' * len(test_structure), name='d1')
strand = nupack.TargetStrand([domain], name='s1')
complex_target = nupack.TargetComplex([strand], test_structure, name='c1')
tube = nupack.TargetTube(on_targets={complex_target: 1e-6}, name='t1')
design = nupack.Design(tubes=[tube], model=config)
results = design.run(trials=1)

print('--- NUPACK ---')
print('Results type:', type(results))
if results:
    trial = results[0]
    print('Trial type:', type(trial))
    print('Trial dir:', [d for d in dir(trial) if not d.startswith('_')])
    try:
        print('Trial.to_analysis() result type:', type(trial.to_analysis()))
    except Exception as e:
        print('Trial.to_analysis() failed:', e)
    
    # Try alternate access
    if hasattr(trial, 'strands'):
        print('Trial.strands type:', type(trial.strands))
    if hasattr(trial, 'tubes'):
        print('Trial.tubes type:', type(trial.tubes))

