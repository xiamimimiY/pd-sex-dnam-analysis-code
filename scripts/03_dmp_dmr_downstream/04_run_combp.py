import subprocess

exploratory_command = 'python -m cpv.pipeline -c 4 --dist 500 --seed 1.0e-4 --anno hg19 -p out forcombp.bed'
final_command = 'python -m cpv.pipeline -c 4 --dist 750 --seed 0.05 --anno hg19 -p out1 combp.bed'

subprocess.call(exploratory_command, shell=True)
subprocess.call(final_command, shell=True)

