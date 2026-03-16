import nupack
import sys

def inspect_obj(name, obj):
    print(f'--- {name} ---')
    print('Type:', type(obj))
    print('Dir:', [d for d in dir(obj) if not d.startswith('_')])
    try:
        print('Repr:', repr(obj))
    except:
        pass

# NUPACK
test_structure = '((((....))))'
config = nupack.Model(material='rna')
domain = nupack.Domain('N' * len(test_structure), name='d1')
strand = nupack.TargetStrand([domain], name='s1')
complex_target = nupack.TargetComplex([strand], test_structure, name='c1')
tube = nupack.TargetTube(on_targets={complex_target: 1e-6}, name='t1')
design = nupack.Design(tubes=[tube], model=config)
results = design.run(trials=1)
inspect_obj('NUPACK Trial', results[0])

# LEARNA/Tensorforce - we need to look at the model inside the agent
# Since we are in the runner, we can't easily import the weights here 
# without the full context, but we can look at the model class if we have an agent.

