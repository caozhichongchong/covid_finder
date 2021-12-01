import pandas as pd
import numpy as np
from datetime import datetime
# load dat
df = pd.read_csv('/scratch/users/anniz44/genomes/covid/clinical/allmut.txt',sep='\t')
df.head()
# get all unique variant
print('finished loading df')
df['variant']=['%s\t%s'%(a,b) for a,b in zip(list(df['POS']), list(df['ALT']))]
print('finished pending variant')

# count abundance
variant_abu = []
allvariantlist = list(df['variant'])
f1 = open('/scratch/users/anniz44/genomes/covid/clinical/allmut.variant.abusum.txt','w')
f1.write('POS\tALT\tabu\n')
i = 0
allvariantdict = dict()
for variant in allvariantlist:
    i += 1
    if i % 10000 == 0:
        allvariantdict[variant] = allvariantdict.get(variant,0) + 1
        print(datetime.now(), 'processed %s variants' % (i))

print('finished counting variant')

for variant in allvariantdict:
    variant_abu.append('%s\t%s\n'%(variant,allvariantdict[variant]))

f1.write(''.join(variant_abu))
f1.close()
print('finished outputing counting variant')

# output to specific state
allvariant = set(list(df['state']))
for state in allvariant:
    variant_abu = df[df['state']==state]
    variant_abu.to_csv('/scratch/users/anniz44/genomes/covid/clinical/%smut.txt'%(state),sep='\t')

print('finished outputing state')
