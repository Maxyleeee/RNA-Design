import nupack
import sys

def test():
    try:
        print("NUPACK version:", nupack.__version__)
        structure = "((....))"
        config = nupack.Model(material='rna')
        target = nupack.TargetStrand(structure, name='t1')
        complex_target = nupack.TargetComplex([target], name='c1')
        design = nupack.Design(complexes=[complex_target], model=config)
        
        print("Running design...")
        results = design.run(trials=1)
        if results:
            seq = results[0].to_analysis().strands[0].sequence
            print("SUCCESS:", seq)
        else:
            print("FAILURE: No results returned")
    except Exception as e:
        print("ERROR:", str(e))
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test()
