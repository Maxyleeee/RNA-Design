import nupack
import sys

def test_nupack_design(structure):
    print(f"Testing NUPACK design for: {structure}")
    try:
        # NUPACK 4 API
        config = nupack.Model(material='rna', ensemble='vienna1999')
        # Define the target structure
        target = nupack.TargetStrand(structure, name='t1')
        # Define the complex
        complex_target = nupack.TargetComplex([target], name='c1')
        # Define the design
        design = nupack.Design(complexes=[complex_target], model=config)
        # Run design
        results = design.run(trials=1)
        # Extract sequence
        if results:
            seq = results[0].to_analysis().strands[0].sequence
            print(f"Success! Sequence: {seq}")
            return seq
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == "__main__":
    test_nupack_design("((....))")
