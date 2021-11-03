import glob
import os
input_script = '/scratch/users/anniz44/scripts/covid/SNPcalling/'

def time_sec(realtime):
    timeinsec = 0
    if 'h' in realtime:
        timeinsec += float(realtime.split('h')[0])*360
        realtime = realtime.split('h')[1]
    if 'm' in realtime:
        timeinsec += float(realtime.split('m')[0])*60
        realtime = realtime.split('m')[1]
    if 's' in realtime:
        timeinsec += float(realtime.split('s')[0])
    return timeinsec

alltime = []
allcov = []
for errfile in glob.glob('%s/*.err'%(input_script)):
    filename = os.path.split(errfile)[-1].split('.mapper1.vcf')[0]
    outfile = errfile.replace('.err','.out')
    time_file = []
    cov_file = []
    for lines in open(errfile,'r'):
        if lines.startswith('real'):
            realtime = lines.split('\t')[1].split('\n')[0]
            time_file.append(time_sec(realtime))
        if 'overall alignment rate' in lines:
            cov_file.append(lines.split('% ')[0])
    for lines in open(outfile,'r'):
        if 'queries matched' in lines:
            cov_file.append(
                '%.2f'%(100*float(lines.split(': ')[1].split('/')[0])/float(lines.split('/')[1].split('\n')[0]))
                            )
    if time_file[0]<time_file[1]:
        print('time ',errfile)
    if float(cov_file[0])>float(cov_file[1]):
        print('cov ',errfile)
    alltime.append('%s\t%.2f\t%.2f\n'%(filename,time_file[0],time_file[1]))
    allcov.append('%s\t%s\t%s\n' % (filename, cov_file[0], cov_file[1]))


f1 = open('%s/../alltimesumcovid.txt'%(input_script),'w')
f1.write('sample\tbowtie\tmapper\n')
f1.write(''.join(alltime))
f1.close()

f1 = open('%s/../allcovsumcovid.txt'%(input_script),'w')
f1.write('sample\tbowtie\tmapper\n')
f1.write(''.join(allcov))
f1.close()