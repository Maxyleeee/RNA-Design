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
print('Trial stats:', trial.stats())
print('Trial domains:', trial.domains)
for d, seq in trial.domains.items():
    print(f'Domain {d.name}: {seq}')

try:
    analysis = trial.to_analysis()
    print('Analysis strands:', analysis.strands)
except Exception as e:
    print('to_analysis error:', e)
