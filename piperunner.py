import os
import re
import argparse
import subprocess
parser = argparse.ArgumentParser(description='Supporting script for StructAlign -> LCore pipe',
                                 epilog = 'Usage example: piperunner.py -lp ../LifeCore/ -i text.txt -o pipe_test.txt -f pipe_test.score -g pipe_test.gk')
parser.add_argument('-lp', type=str, help='Path to LifeCore_main.py programme', dest = 'lcore_path')
parser.add_argument('-i', type=str, help='Input .txt file with PDB IDs and positions', dest = 'inpname')
parser.add_argument('-o', type=str, help='Output .txt file for combined PDB structures', dest = 'outname')
parser.add_argument('-f', type=str, help='Output .txt file for score', dest = 'scorename')
parser.add_argument('-g', type=str, help='Output .txt file for alignments', dest = 'gkname')
args = parser.parse_args()

#run Structalign in current directory
structalign_proc = subprocess.Popen([f'{os.getcwd()}/align', '-i', f'{args.inpname}', '-o', f'{args.outname}',
                        '-f', f'{args.scorename}', '-gk', f'{args.gkname}'], stdout=subprocess.PIPE)

#creating 2 lists, first one – for renamed chains from structalign, second – for original chain names
renaming_old = []
renaming_new = []

#open logfile
with open('pipe.log', 'w') as logfile:
    while True:
      #read structalign stdout and write it to log file
      line = structalign_proc.stdout.readline()
      logfile.write(line.decode('utf8') + '\n')
      if not line:
        break
      line = line.strip().decode('utf8').split(' ')
      #Get last lines of StructAlign output and save chain names
      if line[0] == 'Chain':
          renaming_old.append(line[1])
          renaming_new.append(line[7])

#Read StructAlign gkfile and create temporary file with original naming
with open(f'{args.gkname}', 'r') as gkfile, open('gk_renamed.temp', 'w') as otp:
    for line in gkfile:
        line = line.strip()
        if line[0] != '&':
            line = re.split(r'(:+|\.+|\(+|\)+)', line)
            letter = line[6]
            #rename chains
            for i, j in zip(renaming_old, renaming_new):
                if j == letter:
                    line[6] = i
        otp.write(''.join(line) + '\n')

lcore_proc = subprocess.Popen(['python3', f'{args.lcore_path}LifeCore_main.py',
                                '-i', 'gk_renamed.temp', '-p', '-f'], stdout=subprocess.PIPE)

#Open logfile again and add LCore output
with open('pipe.log', 'a') as logfile:
    while True:
      line = lcore_proc.stdout.readline()
      logfile.write(line.decode('utf8') + '\n')
      if not line:
        break

os.remove("gk_renamed.temp")
print('Output of pipe is pipe.log, gk_renamed_LCore.out, gk_renamed.out.fasta, gk_renamed.pdb.out and gk_renamed.pdb.out.spt files')
