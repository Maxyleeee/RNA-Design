
import sys
import os

log_file = '/home/maxyle/RNA/verify_simple_output_v4.txt'
if not os.path.exists(log_file):
    print(f"Log file {log_file} not found.")
    sys.exit(1)

with open(log_file, 'r') as f:
    content = f.read()

def extract_section(start_marker, end_marker=None):
    if start_marker in content:
        start_idx = content.find(start_marker) + len(start_marker)
        if end_marker:
            end_idx = content.find(end_marker, start_idx)
            if end_idx == -1: return content[start_idx:].strip()
            return content[start_idx:end_idx].strip()
        else:
            return content[start_idx:].strip()
    return 'Marker not found'

print('=== LEARNA STDOUT ===')
learna_out = extract_section('DEBUG LEARNA STDOUT:', 'DEBUG LEARNA STDERR:')
for line in learna_out.split('\n'):
    if '|' in line or 'Prediction' in line:
        print(line)

print('\n=== LEARNA STDERR ===')
learna_err = extract_section('DEBUG LEARNA STDERR:', '--- Verifying NUPACK ---')
print(learna_err)

print('\n=== DesiRNA STDOUT ===')
desirna_out = extract_section('DEBUG DesiRNA STDOUT:', 'DEBUG DesiRNA STDERR:')
for line in desirna_out.split('\n'):
    line = line.strip()
    if 'ETA' not in line and 'Designed' not in line and line and len(line) == 12 and set(line).issubset({'A','C','G','U'}):
        print('FOUND SEQ:', line)

print('\n=== DesiRNA STDERR ===')
desirna_err = extract_section('DEBUG DesiRNA STDERR:', '--- Verifying LEARNA ---')
print(desirna_err)
