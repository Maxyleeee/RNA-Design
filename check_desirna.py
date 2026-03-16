import subprocess
import tempfile
import time
import os

structure = '((.(((((((...)))(((...)))..((.((((....))))))))))))'
desirna_dir = '/home/maxyle/RNA/RNA_design_tools/DesiRNA/DesiRNA-main'
desirna_python = desirna_dir + '/venv/bin/python3'
desirna_script = desirna_dir + '/DesiRNA.py'

seq_restr = 'N' * len(structure)
input_content = f'>name\nBenchmark\n>seq_restr\n{seq_restr}\n>sec_struct\n{structure}\n'

with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, dir='/tmp') as tf:
    tf.write(input_content)
    tmp_path = tf.name

cmd = [desirna_python, desirna_script, '-f', tmp_path, '-R', '1', '-e', '10', '-t', '10']
print('Running:', ' '.join(cmd))
res = subprocess.run(cmd, capture_output=True, text=True)
print('Return code:', res.returncode)
print('STDOUT:')
print(res.stdout)
print('STDERR:')
print(res.stderr)

os.unlink(tmp_path)
