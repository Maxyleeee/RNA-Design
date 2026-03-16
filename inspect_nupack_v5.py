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
print('domains:', trial.domains)
for d in trial.domains:
    print('domain:', d, type(d), dir(d))
    print('domain seq:', getattr(d, 'sequence', 'no sequence'))

try:
    analysis = trial.to_analysis()
    print('analysis:', type(analysis), dir(analysis))
    print('analysis strands:', analysis.strands)
except Exception as e:
    print('to_analysis error:', e)
