# Sort lines of an observational dataset in alphabetical order
import pyneb as pn

# Read data
obs_data = 'pne.dat'
obs = pn.Observation(obs_data, fileFormat='lines_in_rows', corrected=True)

out_file=open('out.dat', 'w')
    
# getSortedLines returns lines in sorted order
for line in obs.getSortedLines():
    row = line.label
    for item in line.corrIntens:
        row += ' {0:10.3e}'.format(item)
    out_file.write(row + '\n')

out_file.close()
