#!/usr/bin/env python
#__________________________________________________________________DNA BARCODE DATA ANALYSIS SCRIPT________________________________________________________________
#
#followup script of the Barcode_Script_1.4 to analyze the data of the output txt files.
#it attributes the filenames to the sample, performs 3 normalization steps and the averaging between the different mice
#written by Sabrina Weis, MSc Student of Molceular Biotechnology, internship March 2017
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

import os
from math import sqrt
import time
import operator

print "\n"*5+"====== Script is running ======\n\n"

def standard_deviation(lst):
    """returns the standard deviation of lst"""
    mn = sum(lst)/len(lst)
    v = sum([(e-mn)**2 for e in lst]) / len(lst)
    return sqrt(v)

execfile("C:/Users/Sabrina/Desktop/SabrinaWeis/Sabrina_script/Barcode_Script.conf")

#open the textfile containing the information which filename is attributed to which sample and store it in a dictionary called sample
with open (filename_sample_file) as temp:
    sample=dict(line.strip().split() for line in temp if line.strip())

#open the textifle containing the noramlization values for each varaint depending on the initial mix; values are stored in a dictionary called variant_normalization
with open(variant_normalization_file) as temp:
    variant_normalization=dict(line.strip().split() for line in temp if line.strip())

#if the values are written with a comma, each comma is replaced by a point for each varaint, cause python can only handle points
for k,v in variant_normalization.items():
    variant_normalization[k]=v.replace(",",".",1)

#the values of the variant_normalization dictionary are converted from a string into a float number to enable calculating
variant_normalization={k:float(v) for k, v in variant_normalization.iteritems()}

#the the variant-file is opend and stored in a dictionary called variants
with open(variants_barcode_file) as temp:
    variants=dict(line.strip().split() for line in temp if line.strip())

values=variants.values()

with open(organ_normalization_file) as temp:
    organ_normalization=dict(line.strip().split() for line in temp if line.strip())

for k,v in organ_normalization.items():
    organ_normalization[k]=v.replace(",",".",1)

organ_normalization={k:float(v) for k, v in organ_normalization.iteritems()}

#all filenames of the given folder are extracted and those which ends with txt are stored in the list inputfiles
inputfiles=[]
objects=os.listdir(my_dir)
for i in range(0,len(objects)):
    filename=objects[i]
    if filename.endswith('txt'):
        if filename in sample.keys():
            inputfiles.append(filename)

#the variable results is defined as a dictionary containing a dictionary for each file stored in the list filenames
results={}
for l in range(0,len(inputfiles)):
    results[inputfiles[l]]={}

#the files in the list inputfiles are opened in succession and each line containing a variant stored in variants are storead in the dictionary results
for n in range(0,len(inputfiles)):
    result_file=open(my_dir+inputfiles[n], 'r')
    for i in result_file:
        line=dict(i.strip().split() for v in values if i.startswith(v))
        if len(line) != 0:
            results[inputfiles[n]].update(line)

#the readcount of the variants is converted from a string into a float number to enabale calculations
for n in results.keys():
    results[n]={k:float(v) for k, v in results[n].iteritems()}

#each filename that occurs in results is replaced by the corresponding sample-name defined in the samplefile
for i in results.keys():
    results[sample[i]]=results[i]
    del results[i]

#the last part of the sample-name is extracted that defines the organ where the sample is coming from; all found organs are stored in a list calld sample_lastpart
sample_lastpart=[]
for n in results.keys():
    sample_lastpart.append(n[8:])

sample_lastpart=sorted(set(sample_lastpart), key=sample_lastpart.index)

#a new dictionary called organs is defined, where the the variants with the readscount are attributed to a sample and each sample is attributed to an organ: {organ:{sample:{variants:readcount}}}

organs_sorted={}
for n in sample_lastpart:
    organs_sorted[n]={}
    for i in results.keys():
        parts=i.split("_")
        if n == parts[2]:
            organs_sorted[n].update({i:results[i]})


cDNA={}
gDNA={}

for n in results:
    if 'cDNA' in n:
        cDNA.update({n[5:]:results[n]})
    elif 'gDNA' in n:
        gDNA.update({n[5:]:results[n]})
    else:
        continue

organs_sorted_cDNA={}
for n in sample_lastpart:
    if len(cDNA) != 0:
        organs_sorted_cDNA[n]={}
        for i in cDNA.keys():
            parts=i.split("_")
            if n == parts[1]:
                organs_sorted_cDNA[n].update({i[:2]:cDNA[i]})
    else:
        break

d1=[]
for i in organs_sorted_cDNA:
    if len(organs_sorted_cDNA[i]) == 0:
        d1.append(i)

for n in d1:
    del organs_sorted_cDNA[n]

organs_sorted_gDNA={}
for n in sample_lastpart:
    if len(gDNA) != 0:
        organs_sorted_gDNA[n]={}
        for i in gDNA.keys():
            parts=i.split("_")
            if n == parts[1]:
                organs_sorted_gDNA[n].update({i[:2]:gDNA[i]})
    else:
        break

d2=[]
for i in organs_sorted_gDNA:
    if len(organs_sorted_gDNA[i]) == 0:
        d2.append(i)

for n in d2:
    del organs_sorted_gDNA[n]

output={'cDNA':organs_sorted_cDNA, 'gDNA':organs_sorted_gDNA}

output_averaged={}
output_stddev={}
for n in output:
    output_averaged[n]={}
    output_stddev[n]={}
    for i in output[n]:
        output_averaged[n][i]={}
        output_stddev[n][i]={}
        for v in output[n][i]['M1']:
            output_averaged[n][i][v]={}
            output_stddev[n][i][v]={}
            counts=[]
            mice=output[n][i].keys()
            for x in range(0,len(mice)):
                key=mice[x]
                if v in output[n][i][key]:
                    counts.append(output[n][i][key][v])
                else:
                    continue
            output_averaged_counts=sum(counts)/len(counts)
            output_stddev_counts=standard_deviation(counts)
            output_averaged[n][i][v]=output_averaged_counts
            output_stddev[n][i][v]=output_stddev_counts

output_sum_counts={}
for n in output:
    output_sum_counts[n]={}
    for i in output[n]:
        output_sum_counts[n][i]={}
        for m in output[n][i]:
            output_sum_counts[n][i][m]=sum(output[n][i][m].values())

flowcell_normalized_output={}
for n in output:
    flowcell_normalized_output[n]={}
    for i in output[n]:
        flowcell_normalized_output[n][i]={}
        for m in output[n][i]:
            flowcell_normalized_output[n][i][m]={}
            for v in output[n][i][m]:
                flowcell_normalized_output[n][i][m][v]=output[n][i][m][v]/output_sum_counts[n][i][m]


organ_normalization_firstpart=[]
for n in organ_normalization.keys():
    organ_normalization_firstpart.append(n[0:2])

organ_normalization_firstpart=sorted(set(organ_normalization_firstpart), key=organ_normalization_firstpart.index)

organ_normalization_struc={}
for n in organ_normalization_firstpart:
    organ_normalization_struc[n]={}
    for i in organ_normalization:
        parts=i.split("_")
        if n == parts[0]:
            organ_normalization_struc[n].update({i[len(n)+1:]:organ_normalization[i]})


organ_normalized_output={}
for n in flowcell_normalized_output:
    organ_normalized_output[n]={}
    for i in flowcell_normalized_output[n]:
        organ_normalized_output[n][i]={}
        for m in flowcell_normalized_output[n][i]:
            organ_normalized_output[n][i][m]={}
            for v in flowcell_normalized_output[n][i][m]:
                organ_normalized_output[n][i][m][v]=flowcell_normalized_output[n][i][m][v]*organ_normalization_struc[m][i]

variant_normalized_output={}
for n in organ_normalized_output:
    variant_normalized_output[n]={}
    for i in organ_normalized_output[n]:
        variant_normalized_output[n][i]={}
        for m in organ_normalized_output[n][i]:
            variant_normalized_output[n][i][m]={}
            for v,c in organ_normalized_output[n][i][m].items():
                for y,z in variant_normalization.items():
                    if y ==v:
                        variant_normalized_output[n][i][m][v]=c/z


stddev={}
average={}
for n in variant_normalized_output:
    average[n]={}
    stddev[n]={}
    for i in variant_normalized_output[n]:
        average[n][i]={}
        stddev[n][i]={}
        for v in variant_normalized_output[n][i]['M1']:
            average[n][i][v]={}
            stddev[n][i][v]={}
            counts=[]
            mice=variant_normalized_output[n][i].keys()
            for x in range(0,len(mice)):
                key=mice[x]
                if v in variant_normalized_output[n][i][key]:
                    counts.append(variant_normalized_output[n][i][key][v])
                else:
                    continue
            averaged_counts=sum(counts)/len(counts)
            stddev_counts=standard_deviation(counts)
            average[n][i][v]=averaged_counts
            stddev[n][i][v]=stddev_counts

variant_normalized_output_sum_variants={}
for n in variant_normalized_output:
    variant_normalized_output_sum_variants[n]={}
    for i in variant_normalized_output[n]:
        variant_normalized_output_sum_variants[n][i]={}
        for m in variant_normalized_output[n][i]:
            variant_normalized_output_sum_variants[n][i][m]=sum(variant_normalized_output[n][i][m].values())

variant_comparison={}
for n in variant_normalized_output:
    variant_comparison[n]={}
    for i in variant_normalized_output[n]:
        variant_comparison[n][i]={}
        for m in variant_normalized_output[n][i]:
            variant_comparison[n][i][m]={}
            for v in variant_normalized_output[n][i][m]:
                variant_comparison[n][i][m][v]=variant_normalized_output[n][i][m][v]/variant_normalized_output_sum_variants[n][i][m]

variant_comparison_average={}
variant_comparison_stddev={}
for n in variant_comparison:
    variant_comparison_average[n]={}
    variant_comparison_stddev[n]={}
    for i in variant_comparison[n]:
        variant_comparison_average[n][i]={}
        variant_comparison_stddev[n][i]={}
        for v in variant_comparison[n][i]['M1']:
            variant_comparison_average[n][i][v]={}
            variant_comparison_stddev[n][i][v]={}
            counts=[]
            mice=variant_comparison[n][i].keys()
            for x in range(0,len(mice)):
                key=mice[x]
                if v in variant_comparison[n][i][key]:
                    counts.append(variant_comparison[n][i][key][v])
                else:
                    continue
            averaged_counts=sum(counts)/len(counts)
            stddev_counts=standard_deviation(counts)
            variant_comparison_average[n][i][v]=averaged_counts
            variant_comparison_stddev[n][i][v]=stddev_counts

variant_normalized_output_sum_organs={}
for n in variant_normalized_output:
    variant_normalized_output_sum_organs[n]={}
    if len(variant_normalized_output[n]) == 0:
        continue
    organs=variant_normalized_output[n].keys()
    for m in variant_normalized_output[n][organs[0]]:
        variant_normalized_output_sum_organs[n][m]={}
        for v in variant_normalized_output[n][organs[0]][m]:
            variant_normalized_output_sum_organs[n][m][v]={}
            counts_organ=[]
            for x in range(0,len(organs)):
                if v in variant_normalized_output[n][organs[x]][m]:
                    counts_organ.append(variant_normalized_output[n][organs[x]][m][v])
                else:
                    continue
            variant_normalized_output_sum_organs[n][m][v]=sum(counts_organ)
            if variant_normalized_output_sum_organs[n][m][v] == 0:
                variant_normalized_output_sum_organs[n][m][v]=1

organ_comparison={}
for n in variant_normalized_output:
    organ_comparison[n]={}
    for i in variant_normalized_output[n]:
        organ_comparison[n][i]={}
        for m in variant_normalized_output[n][i]:
            organ_comparison[n][i][m]={}
            for v in variant_normalized_output[n][i][m]:
                organ_comparison[n][i][m][v]=variant_normalized_output[n][i][m][v]/variant_normalized_output_sum_organs[n][m][v]


organ_comparison_average={}
organ_comparison_stddev={}
for n in organ_comparison:
    organ_comparison_average[n]={}
    organ_comparison_stddev[n]={}
    for i in organ_comparison[n]:
        organ_comparison_average[n][i]={}
        organ_comparison_stddev[n][i]={}
        for v in organ_comparison[n][i]['M1']:
            organ_comparison_average[n][i][v]={}
            organ_comparison_stddev[n][i][v]={}
            counts=[]
            mice=organ_comparison[n][i].keys()
            for x in range(0,len(mice)):
                key=mice[x]
                if v in organ_comparison[n][i][key]:
                    counts.append(organ_comparison[n][i][key][v])
                else:
                    continue
            averaged_counts=sum(counts)/len(counts)
            stddev_counts=standard_deviation(counts)
            organ_comparison_average[n][i][v]=averaged_counts
            organ_comparison_stddev[n][i][v]=stddev_counts

organ_comparison_struc={}
for n in organ_comparison_average:
    organ_comparison_struc[n]={}
    if len(organ_comparison_average[n]) == 0:
        continue
    for v in range(0,len(organ_comparison_average[n][organ_comparison_average[n].keys()[0]])):
        organ_comparison_struc[n][organ_comparison_average[n][organ_comparison_average[n].keys()[0]].keys()[v]]={}
        for i in organ_comparison_average[n]:
            organ_comparison_struc[n][organ_comparison_average[n][i].keys()[v]][i]=organ_comparison_average[n][i][organ_comparison_average[n][i].keys()[v]]
        v+=1


timestr = time.strftime("%Y-%m-%d_%H-%M-%S")
output_file_read_counts=my_dir+'read_counts_'+timestr+'.txt'
output_file_absolute_values=my_dir+'absolute_values_'+timestr+'.txt'
output_file_variant_comparison=my_dir+'variant_comparison_'+timestr+'.txt'
output_file_organ_comparison=my_dir+'organ_comparison_'+timestr+'.txt'

r=open(output_file_read_counts,'w')
r.write("====== Generated with Barcode Script Analysis ======\n")
for i in sorted(sample_lastpart):
    if i in output['cDNA']:
        r.write("\n\n"+i+"_cDNA")
    else:
        r.write("\n\n"+"***")
    if i in output['gDNA']:
        mice_gDNA=len(output['gDNA'][i])
        r.write("\t  "*(mice_gDNA+3)+i+"_gDNA\n")
    else:
        mice_cDNA=len(output['cDNA'][i])
        r.write("\t  "*(mice_cDNA+3)+"***\n")
    if i in output['cDNA']:
        r.write("\nVariants\t")
        for x in sorted(output['cDNA'][i].keys()):
            r.write(str(x)+"\t")
        r.write("Average\tStandard_deviation\t")
    else:
        r.write("\n***\t"+"***\t"*(mice_gDNA+2))
    if i in output['gDNA']:
        r.write("Variants\t")
        for x in sorted(output['gDNA'][i].keys()):
            r.write(str(x)+"\t")
        r.write("Average\tStandard_deviation\t")
    else:
        r.write("***\t"*(mice_cDNA+3))
    if i in output['cDNA']:
        output_averaged_sorted_keys=sorted(output_averaged['cDNA'][i], key=output_averaged['cDNA'][i].get, reverse=True)
        p=0
        for v in output_averaged_sorted_keys:
            r.write("\n"+v+"\t")
            for m in sorted(output['cDNA'][i]):
                r.write(str(output['cDNA'][i][m][v])+"\t")
            r.write(str(output_averaged['cDNA'][i][v])+"\t"+str(output_stddev['cDNA'][i][v])+"\t")
            if i in output['gDNA']:
                gDNA_output_averaged_sorted_keys=sorted(output_averaged['gDNA'][i], key=output_averaged['gDNA'][i].get, reverse=True)
                r.write(gDNA_output_averaged_sorted_keys[p]+"\t")
                for m in sorted(output['gDNA'][i]):
                    r.write(str(output['gDNA'][i][m][gDNA_output_averaged_sorted_keys[p]])+"\t")
                r.write(str(output_averaged['gDNA'][i][gDNA_output_averaged_sorted_keys[p]])+"\t"+str(output_stddev['gDNA'][i][gDNA_output_averaged_sorted_keys[p]])+"\t")
            else:
                r.write("***\t"*(mice_cDNA+3))
            p+=1
    else:
        if i in output['gDNA']:
            gDNA_output_averaged_sorted_keys=sorted(output_averaged['gDNA'][i], key=output_averaged['gDNA'][i].get, reverse=True)
            for g in gDNA_output_averaged_sorted_keys:
                r.write("\n***\t"+"***\t"*(mice_gDNA+2))
                r.write(g+"\t")
                for m in sorted(output['gDNA'][i]):
                    r.write(str(output['gDNA'][i][m][g])+"\t")
                r.write(str(output_averaged['gDNA'][i][g])+"\t"+str(output_stddev['gDNA'][i][g])+"\t")

r.write("\n")
r.close()

av=open(output_file_absolute_values,'w')
av.write("====== Generated with Barcode Script Analysis ======\n\n")
for i in sorted(sample_lastpart):
    if i in sorted(variant_normalized_output['cDNA']):
        av.write("\n\n"+i+"_cDNA")
    else:
        av.write("\n\n"+"***")
    if i in sorted(variant_normalized_output['gDNA']):
        mice_gDNA=len(variant_normalized_output['gDNA'][i])
        av.write("\t  "*(mice_gDNA+3)+i+"_gDNA\n")
    else:
        mice_cDNA=len(variant_normalized_output['cDNA'][i])
        av.write("\t  "*(mice_cDNA+3)+"***\n")
    if i in sorted(variant_normalized_output['cDNA']):
        av.write("\nVariants\t")
        for x in sorted(variant_normalized_output['cDNA'][i].keys()):
            av.write(str(x)+"\t")
        av.write("Average\tStandard_deviation\t")
    else:
        av.write("\n***\t"+"***\t"*(mice_gDNA+2))
    if i in sorted(variant_normalized_output['gDNA']):
        av.write("Variants\t")
        for x in sorted(variant_normalized_output['gDNA'][i].keys()):
            av.write(str(x)+"\t")
        av.write("Average\tStandard_deviation\t")
    else:
        av.write("***\t"*(mice_cDNA+3))
    if i in sorted(variant_normalized_output['cDNA']):
        average_sorted_keys=sorted(average['cDNA'][i], key=average['cDNA'][i].get, reverse=True)
        p=0
        for v in average_sorted_keys:
            av.write("\n"+v+"\t")
            for m in sorted(variant_normalized_output['cDNA'][i]):
                av.write(str(variant_normalized_output['cDNA'][i][m][v])+"\t")
            av.write(str(average['cDNA'][i][v])+"\t"+str(stddev['cDNA'][i][v])+"\t")
            if i in sorted(variant_normalized_output['gDNA']):
                gDNA_average_sorted_keys=sorted(average['gDNA'][i], key=average['gDNA'][i].get, reverse=True)
                av.write(gDNA_average_sorted_keys[p]+"\t")
                for m in sorted(variant_normalized_output['gDNA'][i]):
                    av.write(str(variant_normalized_output['gDNA'][i][m][gDNA_average_sorted_keys[p]])+"\t")
                av.write(str(average['gDNA'][i][gDNA_average_sorted_keys[p]])+"\t"+str(stddev['gDNA'][i][gDNA_average_sorted_keys[p]])+"\t")
            else:
                av.write("***\t"*(mice_cDNA+3))
            p+=1
    else:
        if i in sorted(variant_normalized_output['gDNA']):
            gDNA_average_sorted_keys=sorted(average['gDNA'][i], key=average['gDNA'][i].get, reverse=True)
            for g in gDNA_average_sorted_keys:
                av.write("\n***\t"+"***\t"*(mice_gDNA+2))
                av.write(g+"\t")
                for m in sorted(variant_normalized_output['gDNA'][i]):
                    av.write(str(variant_normalized_output['gDNA'][i][m][g])+"\t")
                av.write(str(average['gDNA'][i][g])+"\t"+str(stddev['gDNA'][i][g])+"\t")

av.write("\n")
av.close()

vc=open(output_file_variant_comparison,'w')
vc.write("====== Generated with Barcode Script Analysis ======\n\n")
for i in sorted(sample_lastpart):
    if i in sorted(variant_comparison['cDNA']):
        vc.write("\n\n"+i+"_cDNA")
    else:
        vc.write("\n\n"+"***")
    if i in sorted(variant_comparison['gDNA']):
        mice_gDNA=len(variant_comparison['gDNA'][i])
        vc.write("\t  "*(mice_gDNA+3)+i+"_gDNA\n")
    else:
        mice_cDNA=len(variant_comparison['cDNA'][i])
        vc.write("\t  "*(mice_cDNA+3)+"***\n")
    if i in sorted(variant_comparison['cDNA']):
        vc.write("\nVariants\t")
        for x in sorted(variant_comparison['cDNA'][i].keys()):
            vc.write(str(x)+"\t")
        vc.write("Average\tStandard_deviation\t")
    else:
        vc.write("\n***\t"+"***\t"*(mice_gDNA+2))
    if i in sorted(variant_comparison['gDNA']):
        vc.write("Variants\t")
        for x in sorted(variant_comparison['gDNA'][i].keys()):
            vc.write(str(x)+"\t")
        vc.write("Average\tStandard_deviation\t")
    else:
        vc.write("***\t"*(mice_cDNA+3))
    if i in sorted(variant_comparison['cDNA']):
        variant_comparison_average_sorted_keys=sorted(variant_comparison_average['cDNA'][i], key=variant_comparison_average['cDNA'][i].get, reverse=True)
        p=0
        for v in variant_comparison_average_sorted_keys:
            vc.write("\n"+v+"\t")
            for m in sorted(variant_comparison['cDNA'][i]):
                vc.write(str(variant_comparison['cDNA'][i][m][v])+"\t")
            vc.write(str(variant_comparison_average['cDNA'][i][v])+"\t"+str(variant_comparison_stddev['cDNA'][i][v])+"\t")
            if i in sorted(variant_comparison['gDNA']):
                gDNA_variant_comparison_average_sorted_keys=sorted(variant_comparison_average['gDNA'][i], key=variant_comparison_average['gDNA'][i].get, reverse=True)
                vc.write(gDNA_variant_comparison_average_sorted_keys[p]+"\t")
                for m in sorted(variant_comparison['gDNA'][i]):
                    vc.write(str(variant_comparison['gDNA'][i][m][gDNA_variant_comparison_average_sorted_keys[p]])+"\t")
                vc.write(str(variant_comparison_average['gDNA'][i][gDNA_variant_comparison_average_sorted_keys[p]])+"\t"+str(variant_comparison_stddev['gDNA'][i][gDNA_variant_comparison_average_sorted_keys[p]])+"\t")
            else:
                vc.write("***\t"*(mice_cDNA+3))
            p+=1
    else:
        if i in sorted(variant_comparison['gDNA']):
            gDNA_variant_comparison_average_sorted_keys=sorted(variant_comparison_average['gDNA'][i], key=variant_comparison_average['gDNA'][i].get, reverse=True)
            for g in gDNA_variant_comparison_average_sorted_keys:
                vc.write("\n***\t"+"***\t"*(mice_gDNA+2))
                vc.write(g+"\t")
                for m in sorted(variant_comparison['gDNA'][i]):
                    vc.write(str(variant_comparison['gDNA'][i][m][g])+"\t")
                vc.write(str(variant_comparison_average['gDNA'][i][g])+"\t"+str(variant_comparison_stddev['gDNA'][i][g])+"\t")

vc.write("\n")
vc.close()


oc=open(output_file_organ_comparison,'w')
oc.write("====== Generated with Barcode Script Analysis ======\n\n")
if len(organ_comparison_struc['cDNA']) != 0:
    key1='cDNA'
else:
    key1='gDNA'
for v in sorted(organ_comparison_struc[key1].keys()):
    if len(organ_comparison['cDNA']) != 0:
        oc.write("\n\n"+v+"_cDNA")
    else:
        oc.write("\n\n***")
    if len(organ_comparison['gDNA']) != 0:
        mice_gDNA=len(organ_comparison['gDNA'][organ_comparison['gDNA'].keys()[0]])
        oc.write("\t"*(mice_gDNA+3)+v+"_gDNA\n")
    else:
        mice_cDNA=len(organ_comparison['cDNA'][organ_comparison['cDNA'].keys()[0]])
        oc.write("\t"*(mice_cDNA+3)+"***\n")
    if len(organ_comparison['cDNA']) != 0:
        oc.write("\nOrgans\t")
        for x in sorted(organ_comparison['cDNA'][organ_comparison['cDNA'].keys()[0]].keys()):
            oc.write(str(x)+"\t")
        oc.write("Average\tStandard_deviation\t")
    else:
        oc.write("\n***\t"+"****\t"*(mice_gDNA+2))
    if len(organ_comparison['gDNA']) != 0:
        oc.write("Organs\t")
        for y in sorted(organ_comparison['gDNA'][organ_comparison['gDNA'].keys()[0]].keys()):
            oc.write(str(y)+"\t")
        oc.write("Average\tStandard_deviation\t")
    else:
        oc.write("****\t"*(mice_cDNA+3))
    if len(organ_comparison['cDNA']) != 0:
        organ_comparison_average_sorted_keys=sorted(organ_comparison_struc['cDNA'][v], key=organ_comparison_struc['cDNA'][v].get, reverse=True)
        length_cDNA=len(organ_comparison_average_sorted_keys)
    else:
        length_cDNA=0
    if len(organ_comparison['gDNA']) != 0:
        gDNA_organ_comparison_average_sorted_keys=sorted(organ_comparison_struc['gDNA'][v], key=organ_comparison_struc['gDNA'][v].get, reverse=True)
        length_gDNA=len(gDNA_organ_comparison_average_sorted_keys)
    else:
        length_gDNA=0
    if length_cDNA >= length_gDNA:
        rounds=length_cDNA
    else:
        rounds=length_gDNA
    for i in range(0,rounds):
        if len(organ_comparison['cDNA']) != 0:
            if len(organ_comparison_average_sorted_keys)-1 >= i:
                oc.write("\n"+organ_comparison_average_sorted_keys[i]+"\t")
                for m in sorted(organ_comparison['cDNA'][organ_comparison_average_sorted_keys[i]]):
                    oc.write(str(organ_comparison['cDNA'][organ_comparison_average_sorted_keys[i]][m][v])+"\t")
                oc.write(str(organ_comparison_average['cDNA'][organ_comparison_average_sorted_keys[i]][v])+"\t"+str(organ_comparison_stddev['cDNA'][organ_comparison_average_sorted_keys[i]][v])+"\t")
                if len(organ_comparison['gDNA'])!= 0:
                    if len(gDNA_organ_comparison_average_sorted_keys)-1 >= i:
                        oc.write(gDNA_organ_comparison_average_sorted_keys[i]+"\t")
                        for m in sorted(organ_comparison['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]]):
                            oc.write(str(organ_comparison['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][m][v])+"\t")
                        oc.write(str(organ_comparison_average['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][v])+"\t"+str(organ_comparison_stddev['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][v])+"\t")
                    else:
                        oc.write("***\t"*(mice_cDNA+3))
                else:
                    oc.write("***\t"*(mice_cDNA+3))
            else:
                if len(organ_comparison['gDNA']) != 0:
                    if len(gDNA_organ_comparison_average_sorted_keys)-1 >= i:
                        oc.write("\n***\t"+"****\t"*(mice_gDNA+2))
                        oc.write(gDNA_organ_comparison_average_sorted_keys[i]+"\t")
                        for m in sorted(organ_comparison['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]]):
                            oc.write(str(organ_comparison['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][m][v])+"\t")
                        oc.write(str(organ_comparison_average['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][v])+"\t"+str(organ_comparison_stddev['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][v])+"\t")
        else:
            if len(gDNA_organ_comparison_average_sorted_keys)-1 >= i:
                oc.write("\n***\t"+"****\t"*(mice_gDNA+2))
                oc.write(gDNA_organ_comparison_average_sorted_keys[i]+"\t")
                for m in sorted(organ_comparison['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]]):
                    oc.write(str(organ_comparison['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][m][v])+"\t")
                oc.write(str(organ_comparison_average['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][v])+"\t"+str(organ_comparison_stddev['gDNA'][gDNA_organ_comparison_average_sorted_keys[i]][v])+"\t")

oc.write("\n")
oc.close()

print "====== Script completed! ======"