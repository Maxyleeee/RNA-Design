import sys
import os

# Add LEARNA to path
sys.path.append(os.path.join(os.getcwd(), 'RNA_design_tools/LEARNA/learna_tools'))

from runner import Runner
import argparse

def test_learna():
    print("Testing LEARNA restore...")
    parser = argparse.ArgumentParser()
    parser.add_argument('--target_structure', type=str, default="((....))")
    parser.add_argument('--restore_path', type=str, default="RNA_design_tools/LEARNA/models/224_0_1")
    parser.add_argument('--num_episodes', type=int, default=1)
    parser.add_argument('--num_steps', type=int, default=10)
    args = parser.parse_args([])
    
    try:
        runner = Runner(args)
        runner.run()
        print("SUCCESS")
    except Exception as e:
        print("ERROR:", str(e))
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_learna()
