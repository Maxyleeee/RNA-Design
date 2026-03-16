import nupack
test_structure = '((((....))))'
config = nupack.Model()
domain = nupack.Domain('N' * len(test_structure), name='d1')
strand = nupack.TargetStrand([domain], name='s1')
complex_target = nupack.TargetComplex([strand], test_structure, name='c1')
tube = nupack.TargetTube(on_targets={complex_target: 1e-6}, name='t1')
dsn = nupack.Design(tubes=[tube], model=config)
results = dsn.run(trials=1)

print('Results type:', type(results))
print('Results dir:', [d for d in dir(results) if not d.startswith('_')])

try:
    analysis = results.to_analysis()
    print('Results.to_analysis() succeeded!')
    print('Analysis type:', type(analysis))
    # print('Analysis strands:', analysis.strands)
except Exception as e:
    print('Results.to_analysis() failed:', e)
