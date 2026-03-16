import nupack
import sys

test_structure = '((((....))))'
config = nupack.Model(material='rna')
domain = nupack.Domain('N' * len(test_structure), name='d1')
strand = nupack.TargetStrand([domain], name='s1')
complex_target = nupack.TargetComplex([strand], test_structure, name='c1')
tube = nupack.TargetTube(on_targets={complex_target: 1e-6}, name='t1')
design = nupack.Design(tubes=[tube], model=config)
results = design.run(trials=1)

trial = results[0]
print('Trial to_analysis type:', type(trial.to_analysis))
try:
    analysis = trial.to_analysis
    print('Analysis result type:', type(analysis))
    print('Analysis dir:', [d for d in dir(analysis) if not d.startswith('_')])
except Exception as e:
    print('Accessing trial.to_analysis failed:', e)

try:
    json_data = trial.to_json()
    print('JSON keys:', json_data.keys() if hasattr(json_data, 'keys') else 'not a dict')
    print('JSON content excerpt:', str(json_data)[:200])
except Exception as e:
    print('trial.to_json() failed:', e)

