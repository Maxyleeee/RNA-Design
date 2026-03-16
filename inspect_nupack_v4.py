import nupack
test_structure = '((((....))))'
config = nupack.Model()
domain = nupack.Domain('N' * len(test_structure), name='d1')
strand = nupack.TargetStrand([domain], name='s1')
complex_target = nupack.TargetComplex([strand], test_structure, name='c1')
tube = nupack.TargetTube(on_targets={complex_target: 1e-6}, name='t1')
dsn = nupack.Design(tubes=[tube], model=config)
results = dsn.run(trials=1)

trial = results[0]
print('Trial type:', type(trial))
print('Trial dir:', [d for d in dir(trial) if not d.startswith('_')])

# Try to find anything with 'seq' or 'strand'
for d in dir(trial):
    if 'seq' in d.lower() or 'strand' in d.lower() or 'tube' in d.lower():
        print(f'Attribute {d}: {getattr(trial, d)}')

try:
    print('Trial repr:', repr(trial))
except:
    pass
