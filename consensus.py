################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################

import os
import sys
import re
import csv
import math
import glob
import shutil
import httplib
import platform
import subprocess

import argparse
import textwrap
import warnings

from math import sqrt

from src.aligner import cdsAlign

from Bio import AlignIO, Entrez, SeqIO
from Bio.Nexus import Nexus
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped, generic_dna
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Align import MultipleSeqAlignment






def to_list(str_obj): return [x for x in str_obj.lstrip("[").rstrip("]").split(", ")]

def reverser(a1, a2): return (a1, a2)[1], (a1, a2)[0]

def file_is_empty(path): return os.stat(path).st_size==0

def ungap(seqObj): return Seq("".join([x for x in str(seqObj) if x != "-"]), generic_dna)

def min_num_index(num_list): return [i for i, x in enumerate(num_list) if x ==min(num_list)][0]

def pick_longest(inputObject): return [rec for rec in inputObject if len(rec.seq) >= max([len(x.seq) for x in inputObject])][0]

def strfind(s, ch): return [i for i, ltr in enumerate(s) if ltr == ch]

def key_max(dictObj): return max(dictObj.iterkeys(), key=lambda k: dictObj[k])

def combine_dict(a, b): return dict(a.items() + b.items() + [(k, a[k] + b[k]) for k in set(b) & set(a)])

def min_stop(s1, s2, s3): return mean([translate(s1).count("*"), translate(s2).count("*"), translate(s3).count("*")])

def pdist(s1, s2): return float(len([x for x, y in zip(s1, s2) if x == y and str(x) not in ["-", "?", "X"] and str(y) not in ["-", "?", "X"]]))/len([x for x, y in zip(s1, s2) if str(x) not in ["-", "?", "X"] and str(y) not in ["-", "?", "X"]])


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def sort_dict_by_key(d): return {k:d[k] for k in natural_sort(sorted(d.keys(), key=lambda x: x[0]))}



def groupy(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def two_to_one(listObj):
    newList = list()
    for l_obj in listObj:
        for data in l_obj:
            newList.append(data)
    return newList


def mean(lst):
    """calculates mean"""
    sum = 0
    for i in range(len(lst)):
        sum += lst[i]
    return (sum / len(lst))


def stddev(lst):

    """calculates standard deviation"""
    sum = 0
    mn = mean(lst)
    for i in range(len(lst)):
        sum += pow((lst[i]-mn),2)

    return sqrt(sum/(len(lst)-1))




def checkTranscript(geneName):
    fileTr = open("checkTranscript.log", "a+")
    transcriptFiles = glob.glob("Sequences/" + geneName + "/*.fas")
    for filename in transcriptFiles:
        handle = open(filename, "rU")
        record = list(SeqIO.parse(handle, "fasta"))
        for rec in record:
            aaSeq1 = translate(rec.seq[:-3])
            if aaSeq1.count("*") > 0:
                s1 = rec.seq[:-3]
                s2 = Seq("N", generic_dna) + rec.seq[:-3]
                s3 = Seq("NN", generic_dna) + rec.seq[:-3]
                min_stop_check = min_stop(s1, s2, s3)
                
                seq_check = [x for x in [s1, s2, s3] if translate(x).count("*") == min_stop_check][0]
                aaSeq2 = translate(seq_check)
                
                fileTr.write("Misframe sequence in file %s transcript ID - %s. Count1 - %s, Count2 - %s\n" %(filename, rec.id, aaSeq1.count("*"), aaSeq2.count("*")))




def finder1(taxa, human, d, part):
    """
        Deals with alternate splicing issue by removing segment of sequence
        that has extreme number of mismatching codons
        
        taxa = taxa sequence
        human = human sequence
        @d = reccomended -> 4
        @part = reccomended -> 20
        """
    
    stringObj_human = [human[i:i+part] for i in range(len(human)) if i <= len(human)-part]
    stringObj_taxa = [taxa[i:i+part] for i in range(len(taxa)) if i <= len(taxa)-part]
    skips = ["-", "?"]
    store_pos = list()
    
    for (i, h_seq), (j, t_seq) in zip(enumerate(stringObj_human), enumerate(stringObj_taxa)):
        
        store_mismatch = list()
        for m,n in zip(h_seq, t_seq):
            if m != n and (m not in skips and n not in skips):
                store_mismatch.append(m)
        
        if len(store_mismatch) <= d and set([x for x in t_seq]) != "-":
            store_pos.append(i)

    store_pos = [x for x in store_pos if stringObj_taxa[x:x+1] != "-"]

    return [x for x in groupy(store_pos)]




def finder2(taxa, human, d, part):
    
    """
        Deals with alternate splicing issue by removing segment of sequence
        that has extreme number of mismatching codons
        
        taxa = taxa sequence
        human = human sequence
        @d = reccomended -> 4
        @part = reccomended -> 20
        """
    
    stringObj_human = [human[i:i+part] for i in range(len(human)) if i <= len(human)-part]
    stringObj_taxa = [taxa[i:i+part] for i in range(len(taxa)) if i <= len(taxa)-part]
    
    skips = ["-", "?"]
    store_pos = list()
    for (i, h_seq), (j, t_seq) in zip(enumerate(stringObj_human), enumerate(stringObj_taxa)):
        store_mismatch = list()
        for m,n in zip(h_seq, t_seq):
            if m != n and (m not in skips and n not in skips):
                store_mismatch.append(m)
    
        try:
            if float(len(store_mismatch))/len([x for x in t_seq if x != "-" and x != "?"]) <= float(d)/part and set([x for x in t_seq]) != set("-"):
                store_pos.append(i)
        except ZeroDivisionError:
            pass
                
    store_pos = [x for x in store_pos if taxa[x:x+1] != "-"]
    
    range_match = list()
    for x in groupy(store_pos):
        range_match.append(range(x[0]-1, x[1]))
        range_match.append(range(x[1]-1, x[1]+20))
    
    range_match=two_to_one(range_match)
    exact_match = [j for (i, h_obj), (j, t_obj) in zip(enumerate(human), enumerate(taxa)) if t_obj == h_obj or j in range_match]

    return [y for y in [x for x in groupy(exact_match)] if y[1]-y[0] > 2]




def transfer_to(output, listObject_h):
    try:
        os.mkdir(output)
    except:
        pass

    folders = ["Sequences/" + x + "/Aligned/" for x in listObject_h]
    for folder in folders:
        try:
            recursive_overwrite(folder, output + "/" + folder.split("/")[1])
        except IOError:
            continue

    print("Files are in %s directory\n" %output)



def rem_short(recordObj):
    ids = set([x.id.split("|")[0] for x in recordObj])
    newRec = list()
    store_fragments = dict()
    for i, id in enumerate(ids):
        store_fragments[i] = [x for x in recordObj if id in x.id]
        newRec.append(pick_longest(store_fragments[i]))
    return newRec


def merge_dict(d1, d2):
    if d1 == {}:
        return d2
    else:
        newDict = dict()
        for key, val in d1.items():
            newDict[key] = (val)
            for inkey, inval in d2.items():
                if key == inkey:
                    for childKey, childVal in inval.items():
                        newDict[key][childKey] = (childVal)
            
                if inkey not in d1.keys():
                    newDict[inkey] = (inval)

        return newDict



def recursive_overwrite(src, dest, ignore=None):
    if os.path.isdir(src):
        if not os.path.isdir(dest):
            os.makedirs(dest)
        files = os.listdir(src)
        if ignore is not None:
            ignored = ignore(src, files)
        else:
            ignored = set()
        for f in files:
            if f not in ignored:
                recursive_overwrite(os.path.join(src, f),
                                    os.path.join(dest, f),
                                    ignore)
    else:
        shutil.copyfile(src, dest)



def exon_names(logfileName, data_dict_human, record_original_tot, geneName):
    exon_file_name = dict()
    file_folders = ["Sequences/" + geneName + "/Aligned/"]
    logData = open(logfileName, 'a+')
    for folders in file_folders:
        frame = data_dict_human[geneName]["frame"]
        if frame == "-":
            if int(record_original_tot[geneName][0][0].split("-")[-1]) < int(record_original_tot[geneName][-1][0].split("-")[-1]):
                record_original_tot[geneName] = record_original_tot[geneName][::-1]
            else:
                pass
        elif frame == "+":
            if int(record_original_tot[geneName][0][0].split("-")[-1]) > int(record_original_tot[geneName][-1][0].split("-")[-1]):
                record_original_tot[geneName] = record_original_tot[geneName][::-1]
            else:
                pass
        
        files = glob.glob(folders + "*.*")
        for i, rec in enumerate(record_original_tot[geneName]):
            for filename in files:
                if filename.split("/")[-1][:-4] in rec[0] or filename.split(":")[-1][:-4] in rec[0]:
                    handle = open(filename, "rU")
                    record = list(SeqIO.parse(handle, "nexus"))
                    handle.close()
                    for j, rec in enumerate(record):
                        if "gi|" in rec.id:
                            record[j].id = "Homo_sapiens_CCDS|CCDS"
                    with open("/".join(filename.split("/")[:-1]) + "/exon" + str(i+1) + ".nex", "w") as fp:
                        SeqIO.write(record, fp, "nexus")
                    os.remove(filename)
                    exon_file_name["exon" + str(i+1)] = (filename)
                    logData.write("%s - exon%s\n"%(filename, i))
                    break

    logData.close()

    return exon_file_name


def tblastnWrapper(geneName, recordX, message):
    print message
    
    store_blasted = dict()
    
    gene_files = dict()
    
    with open("Sequences/" + geneName + "/CCDS.fas", "w") as fp:
        SeqIO.write(recordX, fp, "fasta")
    
    taxa_seq_files = [x for x in glob.glob("Sequences/" + geneName + "/*.fas") if "CCDS" not in x]
    gene_files[geneName] = taxa_seq_files
    toolbar_width = len(taxa_seq_files)
    for fnum, files in enumerate(taxa_seq_files):
        p = str((float(fnum+1)/toolbar_width)*100)[:4]
        sys.stdout.write("\r%s%%" %p)
        sys.stdout.flush()
        
        try:
            os.mkdir(files[:-4])
        except:
            pass

        strore_blast_res = files[:-4] + "/BlastFiles"

        try:
            os.mkdir(strore_blast_res)
        except:
            pass

        blast_to = files[:-4] + "/" + files[:-4].split("/")[-1] + "_blast"
        os.system("makeblastdb -in %s -out %s -dbtype nucl -hash_index" %(files, blast_to))
    
        retData = dict()
        storeHead = dict()

        for i, rec in enumerate(recordX):
            
            retData[rec.id] = dict()
            storeHead[geneName + "_" + str(i)] = rec.id
            
            stringNuc = ''
            for nuc in rec.seq:
                stringNuc = stringNuc + nuc
            
            flagThing = False
            tblastn_cline = NcbitblastnCommandline(cmd='tblastn',
                                              db="%s"%(blast_to),
                                              num_threads=1,
                                              outfmt=5,
                                              seg="no",
                                              out=("%s/Result%s.xml"%(strore_blast_res,i))
                                              )


            try:
                stdout, stderr = tblastn_cline(stdin = stringNuc)
            except:
                continue
                                                 
            """Parse Blast output"""
            with open(strore_blast_res + "/Result" + str(i) + ".xml",'r') as xml:
                cFlag = False; eFlag = False; e_val = list(); negFlag = False; fcount = 0
                hitFlag = False
                for line in xml:
                    if re.search('No hits found', line) == None:
                        """Check if the sequence belong to the group"""
                        cFlag = True
                        
                    if re.search('<Hit_id>', line) != None:
                        hitFlag = True
                                
                    if hitFlag == True:
                        if re.search('<Hit_def>', line) != None:
                            description = line.strip().strip("\n")[9:-10]
                            id = line.strip().strip("\n")[9:-10]
                            retData[rec.id][id] = dict()
                            retData[rec.id][id]["description"] = description
                        if re.search('<Hsp_query-from>', line) != None:
                            retData[rec.id][id]["qstart"] = int(line.strip().rstrip().strip("<Hsp_query-from>").strip('</'))
                        if re.search('<Hsp_query-to>', line) != None:
                            retData[rec.id][id]["qstop"] = int(line.strip().rstrip().strip("<Hsp_query-to>").strip('</'))
                        if re.search('<Hsp_hit-from>', line) != None:
                            retData[rec.id][id]["start"] = int(line.strip().rstrip().strip("<Hsp_hit-from>").strip('</'))
                        if re.search('<Hsp_hit-to>', line) != None:
                            retData[rec.id][id]["stop"] = int(line.strip().rstrip().strip("<Hsp_hit-to>").strip('</'))
                        if re.search('<Hsp_evalue>', line) != None:
                            retData[rec.id][id]["eval"] = float(line.strip().rstrip().strip("<Hsp_evalue>").strip('</'))
                        if re.search('<Hsp_query-frame>', line) != None:
                            retData[rec.id][id]["frame"] = int(line.strip().rstrip().strip("<Hsp_query-frame>").strip('</'))
                            if retData[rec.id][id]["frame"] < 0:
                                retData[rec.id][id]["start"], retData[rec.id][id]["stop"] = reverser(retData[rec.id][id]["start"], retData[rec.id][id]["stop"])
                        if re.search('<Hsp_hseq>', line) != None:
                            retData[rec.id][id]["seq"] = line.strip().rstrip().strip("<Hsp_hseq>").strip('</')
                            hitFlag = False


            store_blasted = merge_dict(store_blasted, retData)

    retData = store_blasted
    
    if not storeHead:
        storeHead = dict()

    return retData, storeHead, gene_files

  
def collect_ccds_record(listObject, data_dict, rev=True):
    orderCCDS = dict()
    record_with_frame = dict()
    record_original = dict()
    new_gene_list = dict()
    for geneName in listObject:
        record_with_frame[geneName] = list()
        record_original[geneName] = list()
        
        try:
            ccds_object = data_dict[geneName]
        except KeyError:
            continue
        
        if rev == True:
            ccds_positions = sorted([(int(x.split("-")[0]), int(x.split("-")[1])) for x in ccds_object["pos"]])[::-1]
        else:
            ccds_positions = sorted([(int(x.split("-")[0]), int(x.split("-")[1])) for x in ccds_object["pos"]])
        
        orderCCDS[geneName] = ccds_positions
        
        remaining = 0
        first_flag = False
        
        for seq_coord in ccds_positions:
            
            flagThing = False

            while flagThing != True:
                try:
                    Entrez.email = "sendambuj@gmail.com"
                    handle_in = Entrez.efetch(db="nucleotide",
                                       id=ccds_object["id"],
                                       rettype="fasta",
                                       strand=+1,
                                       seq_start=seq_coord[0] + 1,
                                       seq_stop=seq_coord[1] + 1)
                                       
                    record = SeqIO.read(handle_in, "fasta")
                    if "gi:" not in record.id:
                        record.id = "gi|" + record.id

                    flagThing = True
        
                except (ValueError, IOError, httplib.HTTPException):
                    print "Something is not right. Program is not able to fetch sequences from NCBI."
                    continue
        
        
            
            if first_flag == False:
                if len(record.seq)%3 != 0:
                    if rev == True:
                        sequenceObj = record.seq.reverse_complement() + Seq("N"*(3 - len(record.seq)%3), generic_dna)
                        record_original[geneName].append([record.id, str(sequenceObj)])
                    else:
                        sequenceObj = record.seq + Seq("N"*(3 - len(record.seq)%3), generic_dna)
                        record_original[geneName].append([record.id, str(sequenceObj)])
                    
                    inseq = translate(sequenceObj)
                    if "*" in inseq:
                        new_gene_list.append(geneName)
                        break
                    
                    remaining = len(record.seq)%3
                    record.seq = inseq
                else:
                    if rev == True:
                        sequenceObj = record.seq.reverse_complement()
                        record_original[geneName].append([record.id, str(sequenceObj)])
                    else:
                        sequenceObj = record.seq
                        record_original[geneName].append([record.id, str(sequenceObj)])
                    
                    remaining = len(record.seq)%3
                    record.seq = translate(sequenceObj)
                        
                first_flag = True
        
            elif first_flag == True:
                if remaining != 0:
                    if rev == True:
                        sequenceObj = Seq("N"*(remaining), generic_dna) + record.seq.reverse_complement()
                    else:
                        sequenceObj = Seq("N"*(remaining), generic_dna) + record.seq
                    
                    remaining = len(sequenceObj)%3
                
                    if len(sequenceObj)%3 != 0:
                        sequenceObj = sequenceObj + Seq("N"*(3-len(sequenceObj)%3), generic_dna)
                    
                    record_original[geneName].append([record.id, str(sequenceObj)])
                    inseq = translate(sequenceObj)
                    record.seq = inseq

                elif len(record.seq)%3 != 0:
                    if rev == True:
                        sequenceObj = record.seq.reverse_complement() + Seq("N"*(3 - len(record.seq)%3), generic_dna)
                    else:
                        sequenceObj = record.seq + Seq("N"*(3 - len(record.seq)%3), generic_dna)
                    
                    record_original[geneName].append([record.id, str(sequenceObj)])
                    inseq = translate(sequenceObj)
                    remaining = len(record.seq)%3
                    record.seq = inseq
                else:
                    if rev == True:
                        sequenceObj = record.seq.reverse_complement()
                    else:
                        sequenceObj = record.seq
                
                    record_original[geneName].append([record.id, str(sequenceObj)])
                    remaining = len(sequenceObj)%3
                    record.seq = translate(sequenceObj)
                    
            if record.seq.count("*") > 1:
                new_gene_list.append(geneName)
                break
            
            record_with_frame[geneName].append(record)
            handle_in.close()


    return record_with_frame, set(new_gene_list), orderCCDS, record_original


def collect_consensus_record(gene_fileList):
    record_consensus = dict()
    gene_fileDict = dict()
    for x in gene_fileList:
        gene_fileDict[x.keys()[0]] = (x.values()[0])
    
    for geneName, files in gene_fileDict.items():
        record_consensus[geneName] = dict()
        for file_obj in files:
            handle = open(file_obj, "r")
            record_consensus_Obj = list(SeqIO.parse(handle, "fasta"))
            record_consensus[geneName][file_obj[:-4].split("/")[-1]] = (record_consensus_Obj)
            handle.close()

    return record_consensus


def pullPositive(data_pos, data_neg):
    retData_blast = dict()
    for ccds_id, val in data_pos.items():
        for consensus_seq_id, inval in val.items():
            if inval["frame"] > 0:
                retData_blast[ccds_id] = (val)

    for ccds_id, val in data_neg.items():
        for consensus_seq_id, inval in val.items():
            if inval["frame"] > 0:
                retData_blast[ccds_id] = (val)

    return retData_blast



def rewriteData(blastRecord, geneName, eval_per_length, consensus_record, record_original_tot=None, record_original_tot_mouse=None):
    
    nContent = dict()

    for i, (key, val) in enumerate(blastRecord.items()):
        recordObj = list()
            
        for rec in record_original_tot[geneName]:
            if key.split("+")[0] in rec[0]:
                ccds_seq = rec[1]
                left_Ns = rec[1][0:3].count("N")
                right_Ns = rec[1][-3:].count("N")
    
    
        nContent[key.split("+")[0].split(":")[1]] = (left_Ns, right_Ns)
        
        if len(rec[1]) > 0 and len(rec[1]) <= 10:
            eval_cut = eval_per_length["0-10"]
        elif len(rec[1]) > 10 and len(rec[1]) <= 20:
            eval_cut = eval_per_length["10-20"]
        else:
            eval_cut = eval_per_length["20-ahead"]

    
        with open("Sequences/" + geneName + "/Results/" + key + ".fas", "w") as fp:
            for inkey, inval in val.items():
                taxa_id = inval["description"]
                
                for taxon, rec_Obj in consensus_record[geneName].items():
                    for rec in rec_Obj:
                        if taxa_id in rec.id and rec.id.split("|")[0] != "Homo_sapiens" and inval["eval"] < eval_cut:
                            if left_Ns == 0 and right_Ns == 0:
                                seqObj = rec.seq[ inval["start"]-1:inval["stop"] ]
                            elif left_Ns != 0 and right_Ns == 0:
                                seqObj = rec.seq[ inval["start"]-1-(3-left_Ns):inval["stop"] ]
                            elif left_Ns == 0 and right_Ns != 0:
                                seqObj = rec.seq[ inval["start"]-1:inval["stop"]+(3-right_Ns) ]
                            else:
                                seqObj = rec.seq[ inval["start"]-1-(3-left_Ns):inval["stop"]+(3-right_Ns) ]
                            
                            recordObj.append(SeqRecord(seqObj, id = taxa_id, description = ""))
        
            recordObj = rem_short(recordObj)
        
            if recordObj != None:
                for rec in record_original_tot[geneName]:
                    if key.split("+")[0].split(":")[1] in rec[0]:
                        recordObj.append(SeqRecord(Seq(rec[1][left_Ns:len(rec[1])-right_Ns], generic_dna), id=rec[0], description = "Homo_sapiens_CCDS"))
                        
                        if record_original_tot_mouse != None:
                            try:
                                mouse_seq = [x.seq for x in record_original_tot_mouse[geneName] if key in x.id][0]
                                recordObj.append(SeqRecord(Seq(mouse_seq[left_Ns:len(rec[1])-right_Ns], generic_dna), id=rec[0], description = "Mus_musculus_CCDS"))
                            except:
                                print("No match found for %s in mouse for %s gene" %(key, geneName))
                        
        
            SeqIO.write(recordObj, fp, "fasta")

    return nContent




def writeData(blastRecord, geneName, consensus_record, eval_per_length, orderCCDS_rev=None, orderCCDS_for=None, record_original_tot=None, record_original_tot_mouse=None):
    nContent = dict()
    for i, (key, val) in enumerate(blastRecord.items()):
        
        if orderCCDS_for != None:
            if geneName in orderCCDS_for.keys():
                for m, object in enumerate(orderCCDS_for[geneName]):
                    if str(object[0]+1) + "-" + str(object[1]+1) in key:
                        key = str(object[0]+1) + "-" + str(object[1]+1)

        else:
            if orderCCDS_rev != None:
                for m, object in enumerate(orderCCDS_rev[geneName]):
                    if str(object[0]+1) + "-" + str(object[1]+1) in key:
                        key = str(object[0]+1) + "-" + str(object[1]+1)
        
        recordObj = list()
        
        try:
            os.mkdir("Sequences/" + geneName + "/Results")
        except:
            pass
        
        for rec in record_original_tot[geneName]:
            if key in rec[0]:
                ccds_seq = rec[1]
                left_Ns = rec[1][0:3].count("N")
                right_Ns = rec[1][-3:].count("N")
    
        if "gi|" in key:
            nContent[key.split(":")[-1]] = (left_Ns, right_Ns)
        else:
            nContent[key] = (left_Ns, right_Ns)

        
        with open("Sequences/" + geneName + "/Results/" + key + ".fas", "w") as fp:
            for inkey, inval in val.items():
                taxa_id = inval["description"]
                
                if len(rec[1]) > 0 and len(rec[1]) <= 10:
                    eval_cut = eval_per_length["0-10"]
                elif len(rec[1]) > 10 and len(rec[1]) <= 20:
                    eval_cut = eval_per_length["10-20"]
                else:
                    eval_cut = eval_per_length["20-ahead"]
                
                for taxon, rec_Obj in consensus_record[geneName].items():
                    for rec in rec_Obj:
                        if taxa_id in rec.id and rec.id.split("|")[0] != "Homo_sapiens" and inval["eval"] < eval_cut:
                            if left_Ns == 0 and right_Ns == 0:
                                seqObj = rec.seq[ inval["start"]-1:inval["stop"] ]
                            elif left_Ns != 0 and right_Ns == 0:
                                seqObj = rec.seq[ inval["start"]-1-(3-left_Ns):inval["stop"] ]
                            elif left_Ns == 0 and right_Ns != 0:
                                seqObj = rec.seq[ inval["start"]-1:inval["stop"]+(3-right_Ns) ]
                            else:
                                seqObj = rec.seq[ inval["start"]-1-(3-left_Ns):inval["stop"]+(3-right_Ns) ]
                            
                            recordObj.append(SeqRecord(seqObj, id = taxa_id, description = ""))

            recordObj = rem_short(recordObj)
            
            if recordObj:
                for rec in record_original_tot[geneName]:
                    if key in rec[0]:
                        recordObj.append(SeqRecord(Seq(rec[1][left_Ns:len(rec[1])-right_Ns], generic_dna), id=rec[0], description = "Homo_sapiens_CCDS"))
                        if record_original_tot_mouse != None:
                            try:
                                mouse_seq = [x.seq for x in record_original_tot_mouse[geneName] if key in x.id][0]
                                recordObj.append(SeqRecord(Seq(mouse_seq[left_Ns:len(rec[1])-right_Ns], generic_dna), id=rec[0], description = "Mus_musculus_CCDS"))
                            except:
                                print("No match found for %s in mouse for %s gene" %(key, geneName))

            
            SeqIO.write(recordObj, fp, "fasta")

    return nContent


def find_empty(listObject):
    emptyGenes = dict()
    for geneName in listObject:
        emptyGenes[geneName] = list()
        files = glob.glob("Sequences/" + geneName + "/Results/*.*")
        for filename in files:
            if file_is_empty(filename) == True:
                emptyGenes[geneName].append(filename)

    return {k: v for k, v in emptyGenes.items() if v != []}


def obj_preObj(ccdsObj, searchObj, record_original_tot_gene):
    collectDict = list()
    nucSeqDict = {x[0]:x[1] for x in record_original_tot_gene}
    stopCodons = ["TAG", "TGA", "TAA"]
    for i, x in enumerate(ccdsObj):
        if searchObj.split("/")[-1] in x.id:
            if nucSeqDict[ccdsObj[0].id][-3:] not in stopCodons:
                if nucSeqDict[x.id][-3:] not in stopCodons:
                    newSeq = x.seq + ccdsObj[i+1].seq
                    idObj = x.id + "+" + ccdsObj[i+1].id + "_lead"

                else:
                    newSeq = ccdsObj[i-1].seq + x.seq
                    idObj = x.id + "+" + ccdsObj[i-1].id + "_lag"
            else:
                if x.seq[-3:] not in stopCodons:
                    newSeq =  x.seq + ccdsObj[i-1].seq
                    idObj = x.id + "+" + ccdsObj[i-1].id + "_lead"
                else:
                    newSeq = ccdsObj[i+1].seq + x.seq
                    idObj = x.id + "+" + ccdsObj[i+1].id + "_lag"
        

            collectDict.append((idObj, newSeq))

    return collectDict


def addExon(ccdsObj, geneName, searchObj, data_dict_human):
    if data_dict_human[geneName]["frame"] == "+":
        for i, x in enumerate(ccdsObj):
            if searchObj.split("/")[-1] in x.id:
                newSeq = x.seq + ccdsObj[i+1].seq
                idObj = x.id + "+" + ccdsObj[i+1].id + "_lead"
    else:
        for i, x in enumerate(ccdsObj):
            if searchObj.split("/")[-1] in x.id:
                newSeq = ccdsObj[i-1].seq + x.seq
                idObj = x.id + "+" + ccdsObj[i-1].id + "_lag"

    collectDict.append((idObj, newSeq))

    return collectDict


def getCCDS_record(data_dict_Obj, listObject):
    data_dict_pos = dict()
    data_dict_neg = dict()
    for key, val in data_dict_Obj.items():
        if val["frame"] == "+":
            data_dict_pos[key] = val
        elif val["frame"] == "-":
            data_dict_neg[key] = val
    
    listObject_pos = [key for key, val in data_dict_pos.items() if key in listObject]
    listObject_neg = [key for key, val in data_dict_neg.items() if key in listObject]
    
    if listObject_pos != []:
        print("\nCollecting records for + frame objects: %s\n" %listObject_pos)
        ccds_record_in, listObject_old_in, orderCCDS_for_in, record_original_for_in = \
                                            collect_ccds_record(listObject, data_dict_pos, rev=False)

    else:
        ccds_record_in = dict()
        orderCCDS_for_in = None
        record_original_for_in = dict()
    
    if listObject_neg != []:
        print("\nCollecting records for - frame objects: %s\n" %listObject_neg)
        ccds_record_neg_in, listObject_new_in, orderCCDS_rev_in, record_original_rev_in = \
                                            collect_ccds_record(listObject_neg, data_dict_Obj)
        
        for obj in ccds_record_neg_in.items():
            ccds_record_in[obj[0]] = (obj[1])
    else:
        orderCCDS_rev_in = None
        orderCCDS_rev_in = dict()
        record_original_rev_in = dict()
    
    try:
        record_original_tot_in = dict(record_original_for_in.items() + record_original_rev_in.items())
    except NameError:
        pass
    
    return orderCCDS_for_in, orderCCDS_rev_in, record_original_tot_in, ccds_record_in




def getCCDS(geneList, CCDSfile):
    fdata = open(CCDSfile, 'r').readlines()
    data_dict = dict()
    for lines in fdata:
        if "Withdrawn" not in lines:
            object = lines.split("\t")
            if object[2] not in data_dict.keys():
                data_dict[object[2]] = dict()
                data_dict[object[2]]["id"] = (object[1])
                data_dict[object[2]]["frame"] = (object[6])
                data_dict[object[2]]["pos"] = (to_list(object[9]))
            elif object[2] in data_dict.keys():
                if len(to_list(object[9])) > len(data_dict[object[2]]["pos"]):
                    data_dict[object[2]]["pos"] = (to_list(object[9]))
                else:
                    continue
    
    return data_dict



##############################################################################################
# Main function
##############################################################################################


def exec_mapping(listObject, tag, match_dict=None):
    listObject_human = listObject
    eval_per_length = {"0-10": 1e-1, "10-20": 1e-3, "20-ahead": 1e-5}
    for geneName in listObject_human:
        for files in glob.glob("Sequences/" + geneName + "/BlastFiles/*.*"):
            if "gi" in files or "Result" in files or "-" in files:
                os.remove(files)
        for files in glob.glob("Sequences/" + geneName + "/Results/*.*"):
            os.remove(files)
        for files in glob.glob("Sequences/" + geneName + "/Aligned/*.*"):
            os.remove(files)

    Entrez.email = "sendambuj@gmail.com"

                
    ##############################################################
    # Human CCDS
    ##############################################################

    if tag == "human":
        data_dict_human = getCCDS(listObject_human, "data/CCDS.txt")
    else:
        data_dict_human = getCCDS(listObject_human, "data/CCDS.current.mouse.txt")
        data_dict_human = {match_dict[key]:val for key, val in data_dict_human.items() if key in listObject_human}
        listObject_human = [match_dict[x] for x in listObject_human]

    orderCCDS_for, orderCCDS_rev, record_original_tot, ccds_record = getCCDS_record(data_dict_human, listObject_human)

    blast_out = dict()
    gene_fileList = list()
                                        
    for geneName in listObject_human:
        message = "Running Blast for " + geneName
        blast_out[geneName], headStorage, gene_files = tblastnWrapper(geneName, ccds_record[geneName], message)
        gene_fileList.append(gene_files)
        print("\n")

    consensus_record = collect_consensus_record(gene_fileList)

    nContentRet = dict()

    for key, blastRecord in blast_out.items():
        nContentRet[key] = writeData(blastRecord, key, consensus_record, eval_per_length, orderCCDS_rev, orderCCDS_for, record_original_tot, record_original_tot_mouse=None)



    ########################################################################################################
    # Rerun tblastn for small CCDS querry showing no hits. Program concatenate neighbour
    # CCDS sequence with the querry sequence to increase the overall querry length.
    # sequence hit region for the newly concatenated CCDS sequence is them removed from
    # the final prealigned sequence file
    ########################################################################################################

    emptyFiles = find_empty(listObject_human)

    if emptyFiles != {}:
        print("Folders with empty elements are: %s" %emptyFiles)

        new_ccds_record = dict()
        for geneName, files in emptyFiles.items():
            new_ccds_record[geneName] = list()
            for filename in files:
                collectEmpty = (obj_preObj(ccds_record[geneName], filename[:-4], record_original_tot[geneName]))
                for objects in collectEmpty:
                    new_ccds_record[geneName].append(SeqRecord(objects[1], id = objects[0]))

        blast_out_new = dict()

        for geneName in emptyFiles.keys():
            message = "Re-running blast for " + geneName
            blast_out_new[geneName], headStorage, gene_files = tblastnWrapper(geneName, new_ccds_record[geneName], message)

        for key, blastRecord in blast_out_new.items():
            nContentRet[key] = combine_dict(nContentRet[key], rewriteData(blastRecord, key, eval_per_length, consensus_record, record_original_tot, record_original_tot_mouse=None))

        for geneName in blast_out_new.keys():
            files = [x for x in glob.glob("Sequences/" + geneName + "/Results/*.*") if "_lead" in x or "_lag" in x]
            prefiles = [x for x in glob.glob("Sequences/" + geneName + "/Results/*.*") if "_lead" not in x and "_lag" not in x]
            for filename in files:
                handle = open(filename, 'rU')
                recordObj_twins = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                if "_lead" in filename:
                    leadElement = filename.split(":")[2][:-9]
                    leadFile = [x for x in prefiles if leadElement in x and "_lead" not in x][0]
                    leadHandle = open(leadFile, 'rU')
                    leadRecord = SeqIO.to_dict(SeqIO.parse(leadHandle, "fasta"))
                    leadRecord = {k.split("|")[0]:v for k, v in leadRecord.items()}
                    leadHandle.close()
                    for i, rec in enumerate(recordObj_twins):
                        if "gi|" in rec.id:
                            continue
                        else:
                            try:
                                s1 = str(leadRecord[rec.id.split("|")[0]].seq[:10]).strip("N")
                                s = str(rec.seq)
                                try:
                                    pos = re.search(str(s1), str(s)).start()
                                    recordObj_twins[i].seq = rec.seq[:pos]
                                except AttributeError:
                                    recordObj_twins[i].seq = rec.seq[:-len(leadRecord[rec.id.split("|")[0]])]
                                    continue
                            except KeyError:
                                continue
        
                elif "_lag" in filename:
                    lagElement = filename.split(":")[2][:-8]
                    lagFile = [x for x in prefiles if lagElement in x and "_lag" not in x][0]
                    lagHandle = open(lagFile, 'rU')
                    lagRecord = SeqIO.to_dict(SeqIO.parse(lagHandle, "fasta"))
                    lagRecord = {k.split("|")[0]:v for k, v in lagRecord.items()}
                    lagHandle.close()
                    for i, rec in enumerate(recordObj_twins):
                        if "gi|" in rec.id:
                            continue
                        else:
                            try:
                                s1 = str(lagRecord[rec.id.split("|")[0]].seq[-10:]).strip("N")
                                s = str(rec.seq)
                                try:
                                    pos = re.search(str(s1), str(s)).end()
                                    recordObj_twins[i].seq = rec.seq[pos:]

                                except AttributeError:
                                    recordObj_twins[i].seq = rec.seq[len(lagRecord[rec.id.split("|")[0]]):]
                                    continue
                            except KeyError:
                                continue


                newRecObject_twins = list()
                for rec in recordObj_twins:
                    if len(rec.seq) != 0:
                        newRecObject_twins.append(rec)
            
                os.remove(filename)

                with open(filename.split("+")[0].split("gi")[0] + filename.split("+")[0].split(":")[-1] + ".fas", "w") as fp:
                    SeqIO.write(newRecObject_twins, fp, "fasta")
                    
                        
            
    return nContentRet, listObject_human, data_dict_human, record_original_tot


##########################################################################################################
# Exon alignments
##########################################################################################################

def align_exon(nContentRet, listObject_h):
    exon_failed = list()
    seqFolders = ["Sequences/" + x + "/" for x in listObject_h]
    logfile = open("logData.log", "a+")
    for folders in seqFolders:
    
        logfile.write("\n\n%s\n" %folders[10:-1])
    
        files = [x for x in glob.glob(folders + "Results/*.fas") if "_lag" not in x and "_lead" not in x]
        for filename in files:
            os.rename(filename, filename.replace("|", "-"))
            filename = filename.replace("|", "-")
        
            try:
                os.mkdir(filename[:strfind(filename, "/")[-2]+1] + "Aligned")
            except:
                pass
        
        
            if "gi-" in filename.split("/")[-1][:-4]:
                fileCheck = filename.split("/")[-1][:-4].split(":")[1]
            else:
                fileCheck = filename.split("/")[-1][:-4]
            
            geneName = folders[10:len(folders)-1]
        
            if fileCheck in nContentRet[geneName].keys():
                left_Ns = nContentRet[geneName][fileCheck][0]*"N"
                right_Ns = nContentRet[geneName][fileCheck][1]*"N"

            handle = open(filename, 'rU')
            record = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            
            checkrec = list()

            with open(filename, "w") as fp:
                print filename
                for i, rec in enumerate(record):
                    if len(rec.seq) == 0:
                        continue
                    
                    record[i].seq = Seq(left_Ns, generic_dna) + rec.seq + Seq(right_Ns, generic_dna)
                    checkrec.append(record[i]) #Remove all the zero length sequences
            
                SeqIO.write(checkrec, fp, "fasta")

            outfile = filename.replace("Results", "Aligned")

            try:
                cdsAlign(filename, outfile, stage="mapper")
            except:
                if "gi" not in filename:
                    exon_failed.append(filename)
                if checkrec != []:
                    with open(outfile[:-3] + "nex", "w") as fp:
                        SeqIO.write(checkrec, fp, "nexus")
                else:
                    continue
            
                continue
            
            handle = open(outfile, 'rU')
            record = list(SeqIO.parse(handle, "fasta"))
            
            handle.close()

            for rec in record:
                if "gi|" in rec.id:
                    positions = [i for i, x in enumerate(rec.seq) if x == "-"]


            # Do not skip insertion
            #positions = []
            #######################
    
            seqRecordObj = list()
            for rec in record:
                sequenceObj = Seq("", generic_dna)
                for i, nuc in enumerate(rec.seq):
                    if i not in positions:
                        sequenceObj = sequenceObj + Seq(str(nuc), generic_dna)

            
                if len([x for x in rec.seq if x != "N" and x != "-"]) != 0:
                    seqRecordObj.append(SeqRecord(Seq(str(sequenceObj), generic_dna), id=rec.id, name=rec.name, description=rec.description))
                    
            ref_seq_record = list()
            human_seq = [x.seq for x in seqRecordObj if "gi|" in x.id][0]
            with open("pdist_log.txt", "a+") as fp:
                fp.write("\n%s\n\n"%folders.split("/")[-1])
                for rec in seqRecordObj:
                    p_distance = pdist(rec.seq, human_seq)
                    if p_distance >= 0.75:
                        ref_seq_record.append(rec)
                    else:
                        pass
            
                    fp.write("%s\t%s\t%s\n" %(filename, rec.id, p_distance))
            
            with open(outfile[:-3] + "nex", "w") as fp:
                SeqIO.write(ref_seq_record, fp, "nexus")

            os.remove(outfile)

    logfile.close()

    return exon_failed


####################################################################################
# Fill gaps by CDS mapping
####################################################################################



def ungap(seqObj): return Seq("".join([x for x in str(seqObj) if x != "-"]), generic_dna)

def seqInitGap(seqObj):
    for i, nuc in enumerate(seqObj):
        if nuc == "-":
            gap_init = i
        else:
            break
    try:
        gap_init
    except UnboundLocalError:
        gap_init = -1
    return gap_init


def findTermGap(fileObj):
    handle = open(fileObj, "rU")
    record = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    for i, nuc in enumerate(record[1].seq):
        if nuc == "-":
            gap_init = i
        else:
            break

    for i, nuc in enumerate(record[1].seq[::-1]):
        if nuc == "-":
            gap_end = len(record[1].seq) - i
        else:
            break

    try:
        gap_init
    except UnboundLocalError:
        gap_init = None

    try:
        gap_end
    except UnboundLocalError:
        gap_end = None

    return gap_init, gap_end, str(record[0].seq)


def trimmer(querry, target, left_Ns, right_Ns, posCheck, init_seq, exon_left_Ns, exon_right_Ns):
    left_Ns = left_Ns + (posCheck%3)*"N"
    querry_start = target.find(querry)
    if querry_start != -1:
        return target[querry_start+len("".join([x for x in querry if x != "-"])):]
    else:
        ss1 = str(init_seq)
        ss2 = str(target)
        
        # Find the position of gapped pos in trimmed sequence
        # Walk left in trimmed sequence by counting number of gaps + number of actual Ns in the gapped sequence
        # Now count the number of remaining nucleotides towards the left in trimmed sequence
        # do 3 - (number of nucs)%3
        
        findPos = ss2.find("".join([x for x in ss1 if x != "-"]))
        if findPos == -1:
            findPos = ss2.find("".join([x for x in ss1 if x != "-"])[:15])

        if findPos == -1:
            writeRecord1 = [SeqRecord(Seq(ss2, generic_dna), id="Seq1"), SeqRecord(Seq("".join([x for x in ss1 if x != "-"]), generic_dna), id="Seq2")]
                
            with open("timmer.fas", "w") as fq:
                SeqIO.write(writeRecord1, fq, "fasta")
                
            if 'Darwin' in platform.system():
                subprocess.call("./src/muscle/muscle -in timmer.fas -out timmer_out.fas", shell=True)
            else:
                subprocess.call("./src/muscle/muscleLinux -in timmer.fas -out timmer_out.fas", shell=True)
                
            handle_fas = open("timmer_out.fas", "rU")
            record_fas = list(SeqIO.parse(handle_fas, "fasta"))
            counter = 0
            for nuc1, nuc2 in zip(record_fas[0].seq, record_fas[1].seq):
                if "-" in nuc2:
                    counter = counter + 1
                else:
                    break

            os.remove("timmer_out.fas")
            os.remove("timmer.fas")
            findPos = counter

        extend_range = findPos - (len(exon_left_Ns) + seqInitGap(ss1) + 1)
        if extend_range < 0:
            extend_range = extend_range + 3

        extNs = 3 - len(ss2[:extend_range])%3
        if extNs == 3:
            extNs = 0

        if extend_range == 0 or findPos == 0:
            extNs = len(exon_left_Ns)

        target = Seq("N"*extNs, generic_dna) + Seq("".join([x for x in str(target) if x != "-"]), generic_dna)
        extNR = 3 - len(target)%3
        targetSeq = SeqRecord(target + Seq(extNR*"N", generic_dna), id="Target")
        querrySeq = SeqRecord(Seq(left_Ns, generic_dna) + Seq(querry, generic_dna) + Seq(right_Ns, generic_dna), id="Querry")
        writeRecord = [targetSeq, querrySeq]
        with open("timmer.fas", "w") as fp:
            SeqIO.write(writeRecord, fp, "fasta")
        
        cdsAlign("timmer.fas", "timmer_out.fas", stage="mapper1")
        gap_init, gap_end, target_ret = findTermGap("timmer_out.fas")
        os.remove("timmer.fas")
        os.remove("timmer_out.fas")
        target = target[extNs:-extNR]
        if gap_end != None:
            return str(target_ret[gap_end-1-len(right_Ns):].lstrip("N").rstrip("N"))
        else:
            return str(target)


def trimmer1(querry, target, left_Ns, right_Ns, posCheck, init_seq, exon_left_Ns, exon_right_Ns):
    left_Ns = left_Ns + (posCheck%3)*"N"
    querry_start = target.find("".join([x for x in querry if x != "-"]))
    if querry_start != -1:
        return target[:querry_start]
    else:
        ss1 = str(init_seq)
        ss2 = str(target)
            
        # Find the position of gapped seq in trimmed sequence
        # Walk left in trimmed sequence by counting number of gaps + number of actual Ns in the gapped sequence
        # Now count the number of remaining nucleotides towards the left in trimmed sequence
        # do 3 - (number of nucs)%3
            
        findPos = ss2.find("".join([x for x in ss1 if x != "-"]))
        if findPos == -1:
            findPos = ss2.find("".join([x for x in ss1 if x != "-"])[:15])
        
        if findPos == -1:
            writeRecord1 = [SeqRecord(Seq(ss2, generic_dna), id="Seq1"), SeqRecord(Seq("".join([x for x in ss1 if x != "-"]), generic_dna), id="Seq2")]
            
            with open("timmer.fas", "w") as fq:
                SeqIO.write(writeRecord1, fq, "fasta")
            
            print "####################################################################"
            
            if 'Darwin' in platform.system():
                subprocess.call("./src/muscle/muscle -in timmer.fas -out timmer_out.fas", shell=True)
            else:
                subprocess.call("./src/muscle/muscleLinux -in timmer.fas -out timmer_out.fas", shell=True)
            
            handle_fas = open("timmer_out.fas", "rU")
            record_fas = list(SeqIO.parse(handle_fas, "fasta"))
            
            os.remove("timmer_out.fas")
            os.remove("timmer.fas")
            
            counter = 0
            
            for nuc1, nuc2 in zip(record_fas[0].seq, record_fas[1].seq):
                if "-" in nuc2:
                    counter = counter + 1
                else:
                    break
            
            findPos = counter
        
        extend_range = findPos - (len(exon_left_Ns) + seqInitGap(ss1) + 1)

        if extend_range < 0:
            extend_range = extend_range + 3

        extNs = 3 - len(ss2[:extend_range])%3
            
        if extNs == 3:
            extNs = 0

        if extend_range == 0 or findPos == 0:
            extNs = len(exon_left_Ns)

        target = Seq("N"*extNs, generic_dna) + Seq("".join([x for x in str(target) if x != "-"]), generic_dna)
        extNR = 3 - len(target)%3
        targetSeq = SeqRecord(target + Seq(extNR*"N", generic_dna), id="Target")
        querrySeq = SeqRecord(Seq(left_Ns, generic_dna) + Seq(querry, generic_dna) + Seq(left_Ns, generic_dna), id="Querry")
        writeRecord = [targetSeq, querrySeq]
        with open("timmer.fas", "w") as fp:
            SeqIO.write(writeRecord, fp, "fasta")
        
        cdsAlign("timmer.fas", "timmer_out.fas", stage="mapper1")
        gap_init, gap_end, target_ret = findTermGap("timmer_out.fas")
        os.remove("timmer.fas")
        os.remove("timmer_out.fas")
        
        target = target[extNs:-extNR]

        if gap_init != None:
            return str(target_ret[:gap_init].lstrip("N").rstrip("N"))
        else:
            return str(target)



def mapper_whole(geneName, nContentRet, ex_file_name, outlier_files):
        
    outlier_files.write("[%s]\n\n" %geneName)
    saveOrig = open("GappedSeq.txt", "a+")
    exonFiles = glob.glob("Sequences/" + geneName + "/Aligned/*.nex")
    exon_dict = dict()
    store_target = list()
    tempStore = dict()
    for i, filename in enumerate(exonFiles):
        
        handle = open(filename, "rU")
        record = list(SeqIO.parse(handle, "nexus"))
        handle.close()
        for rec in record:
            if i == 0:
                exon_dict[str(rec.id).split("|")[0]] = dict()
            try:
                exon_dict[str(rec.id).split("|")[0]][filename.split("/")[-1]] = str(rec.seq)
            except KeyError:
                exon_dict[str(rec.id).split("|")[0]] = dict()
                exon_dict[str(rec.id).split("|")[0]][filename.split("/")[-1]] = str(rec.seq)
            
            if "-" in rec.seq[:3] or "-" in rec.seq[-3:]:
                if filename not in tempStore.keys():
                    tempStore[filename] = dict()
                tempStore[filename][rec.id] = rec.seq
                saveOrig.write("%s -> %s\n" %(rec.id, str(rec.seq)))
                store_target.append([filename, str(rec.id).split("|")[0], str(rec.id).split("|")[1]])

    for files in exonFiles:
        handle = open(files, "rU")
        record = list(SeqIO.parse(handle, "nexus"))
        with open(files[:-4] + ".fas", "w") as fp:
            SeqIO.write(record, fp, "fasta")



    for work_element in store_target:
        
        try:
            transcript_seq_handle = open("Sequences/" + geneName + "/" + work_element[1] + ".fas", "rU")
        except IOError:
            continue
        
        record_transcript = list(SeqIO.parse(transcript_seq_handle, "fasta"))
        for x in record_transcript:
            if work_element[2] in x.id:
                transcript_seq = str(x.seq)

        key_list = natural_sort(exon_dict[work_element[1]].keys())
        flag = False

        # Find number of Ns for exon in work_element[0]

        exon_fileCheck = ex_file_name[work_element[0].split("/")[-1][:-4]].split("/")[-1][:-4]
        if "gi-" in exon_fileCheck:
            exon_fileCheck = exon_fileCheck.split(":")[1]

        if exon_fileCheck in nContentRet[geneName].keys():
            exon_left_Ns = nContentRet[geneName][exon_fileCheck][0]*"N"
            exon_right_Ns = nContentRet[geneName][exon_fileCheck][1]*"N"

        ########################################


        for counter, key in enumerate(key_list):
            fileCheck = ex_file_name[key[:-4]].split("/")[-1][:-4]

            if "gi-" in fileCheck:
                fileCheck = fileCheck.split(":")[1]

            if fileCheck in nContentRet[geneName].keys():
                left_Ns = nContentRet[geneName][fileCheck][0]*"N"
                right_Ns = nContentRet[geneName][fileCheck][1]*"N"
    
            if flag == True:
                val = exon_dict[work_element[1]][key]
                posCheck = 0
                for i, nuc in enumerate(val):
                    if nuc == "-":
                        posCheck = i + 1
                    else:
                        break

                val1 = "".join([x for x in val if x != "-"])
                transcript_seq = trimmer1(val1, transcript_seq, left_Ns, right_Ns, posCheck, exon_dict[work_element[1]][work_element[0].split("/")[-1]], exon_left_Ns, exon_right_Ns)
                break
                    
            if key in work_element[0]:
                flag = True
            
            elif flag == False and "exon" + str(int(key[-5:-4]) + 1) + ".nex" in work_element[0]:
                val = exon_dict[work_element[1]][key]
                fileCheck = ex_file_name[key[:-4]].split("/")[-1][:-4]
                if "gi-" in fileCheck:
                    fileCheck = fileCheck.split(":")[1]

                if fileCheck in nContentRet[geneName].keys():
                    left_Ns = nContentRet[geneName][fileCheck][0]*"N"
                    right_Ns = nContentRet[geneName][fileCheck][1]*"N"
    
                posCheck = 0
                for i, nuc in enumerate(val):
                    if nuc == "-":
                        posCheck = i + 1
                    else:
                        break
                
                val1 = "".join([x for x in val if x != "-"])
                transcript_seq = trimmer(val1, transcript_seq, left_Ns, right_Ns, posCheck, exon_dict[work_element[1]][work_element[0].split("/")[-1]], exon_left_Ns, exon_right_Ns)

        try:
            handle = open(work_element[0][:-4] + ".fas", "rU")
        except IOError:
            continue

        record = list(SeqIO.parse(handle, "fasta"))
        for i, rec in enumerate(record):
            if rec.id.split("|")[0] in work_element[1] and "Homo_sapiens" not in work_element[1]:
                record[i].seq = Seq(transcript_seq, generic_dna)

        with open(work_element[0][:-4] + ".fas", "w") as fp:
            SeqIO.write(record, fp, "fasta")
                
        

    for filename in exonFiles:
        fileCheck = ex_file_name[filename.split("/")[-1][:-4]].split("/")[-1][:-4]
        if "gi-" in fileCheck:
            fileCheck = fileCheck.split(":")[1]

        if fileCheck in nContentRet[geneName].keys():
            left_Ns = nContentRet[geneName][fileCheck][0]*"N"
            right_Ns = nContentRet[geneName][fileCheck][1]*"N"

        try:
            handle = open(filename[:-4] + ".fas", "rU")
        except IOError:
            continue
        
        record = list(SeqIO.parse(handle, "fasta"))
        for rec in record:
            if "Homo_sapiens" in rec.id:
                refSeq = (Seq(left_Ns, generic_dna) + rec.seq + Seq(right_Ns, generic_dna)).translate()

        N_dict = dict()
        with open(filename[:-4] + "withN.fas", "w") as fp:
            
            for i, rec in enumerate(record):
                if "exon1.nex" in filename or filename not in tempStore.keys():
                    record[i].seq = Seq(left_Ns, generic_dna) + Seq("".join([x for x in str(rec.seq) if x != "-"]), generic_dna) + Seq(right_Ns, generic_dna)
                    record[i].seq = ungap(record[i].seq)
                    continue
            
                if rec.id in tempStore[filename].keys():
                    if [filename, str(rec.id).split("|")[0], str(rec.id).split("|")[1]] in store_target:
                        ss1 = str(tempStore[filename][rec.id])
                        ss2 = str(rec.seq)
                        
                        # Find the position of gapped seq in trimmed sequence
                        # Walk left in trimmed sequence by counting number of gaps + number of actual Ns in the gapped sequence
                        # Now count the number of remaining nucleotides towards the left in trimmed sequence
                        # do 3 - (number of nucs)%3
                        
                        findPos = ss2.find("".join([x for x in ss1 if x != "-"]))
                        if findPos == -1:
                            findPos = ss2.find("".join([x for x in ss1 if x != "-"])[:15])
                        
                        if findPos == -1:
                            writeRecord = [record[i], SeqRecord(tempStore[filename][rec.id], id="seq2")]
                            with open("timmer.fas", "w") as fq:
                                SeqIO.write(writeRecord, fq, "fasta")

                            if 'Darwin' in platform.system():
                                subprocess.call("./src/muscle/muscle -in timmer.fas -out timmer_out.fas", shell=True)
                            else:
                                subprocess.call("./src/muscle/muscleLinux -in timmer.fas -out timmer_out.fas", shell=True)
                            
                            handle_fas = open("timmer_out.fas", "rU")
                            record_fas = list(SeqIO.parse(handle_fas, "fasta"))
                            counter = 0
                            try:
                                for nuc in tempStore[filename][rec.id]:
                                    if "-" in nuc:
                                        counter = counter + 1
                                    else:
                                        break
                            except IndexError:
                                if left_Ns == 3:
                                    left_Ns = 0
                                
                                record[i].seq = Seq(left_Ns, generic_dna) + Seq("".join([x for x in str(rec.seq) if x != "-"]), generic_dna) + Seq(right_Ns, generic_dna)
                
                                if record[i].seq[:3].count("N") == 3:
                                    record[i].seq = record[i].seq[3:]
                
                                record[i].seq = ungap(record[i].seq)
                                continue
                            
                            if counter == 0:
                                extNs = len(left_Ns)

                            extNs = 3 - (counter - len(left_Ns))%3
                            if extNs == 3:
                                extNs = 0
                            
                            record[i].seq = Seq(extNs*"N", generic_dna) + Seq("".join([x for x in str(rec.seq) if x != "-"]), generic_dna)
                            extNR = 3 - len(record[i].seq)%3
                            record[i].seq = record[i].seq + Seq(extNR*"N", generic_dna)
                            if record[i].seq[:3].count("N") == 3:
                                record[i].seq = record[i].seq[3:]
                            
                            record[i].seq = ungap(record[i].seq)
                            continue
                    
                        if findPos == seqInitGap(ss1)+1:
                            extNs = len(left_Ns)
                        
                            if extNs == 3:
                                extNs = 0

                            record[i].seq = Seq("N"*extNs, generic_dna) + Seq("".join([x for x in str(rec.seq) if x != "-"]), generic_dna)
                            extNR = 3 - len(record[i].seq)%3
                            record[i].seq = record[i].seq + Seq(extNR*"N", generic_dna)
                            if record[i].seq[:3].count("N") == 3:
                                record[i].seq = record[i].seq[3:]
                            
                            record[i].seq = ungap(record[i].seq)
                            continue

                        extend_range = findPos - (len(left_Ns) + seqInitGap(ss1) + 1)
                        if extend_range < 0:
                            extend_range = extend_range + 300

                        extNs = 3 - len(ss2[:extend_range])%3
                        if extNs == 3:
                            extNs = 0
                        
                        if extend_range == 0:
                            extNs = 0
                        elif findPos == 0:
                            extNs = len(left_Ns)

                        record[i].seq = Seq("N"*extNs, generic_dna) + Seq("".join([x for x in str(rec.seq) if x != "-"]), generic_dna)
                        extNR = 3 - len(record[i].seq)%3
                        record[i].seq = record[i].seq + Seq(extNR*"N", generic_dna)
                        tr_seq = record[i].seq.translate()
                        if tr_seq.count("*") > 1 and (findPos < (len(left_Ns) + seqInitGap(ss1) + 1)):
                            extNs = len(left_Ns)
                            record[i].seq = Seq("N"*extNs, generic_dna) + Seq(ss2, generic_dna)
                            extNR = 3 - len(record[i].seq)%3
                            record[i].seq = record[i].seq + Seq(extNR*"N", generic_dna)

                        translated_seq = record[i].seq[:-6].translate()
                        addN = 0
                        if translated_seq.count("*") > 0:
                            count1 = translated_seq.count("*")
                            s_check1 = Seq("N", generic_dna) + Seq("".join([x for x in str(record[i].seq) if x != "-"]), generic_dna)
                            translated_seq = s_check1.translate()
                            if "*" in translated_seq:
                                count2 = translated_seq.count("*")
                                s_check2 = Seq("N", generic_dna) + Seq("".join([x for x in str(s_check1) if x != "-"]), generic_dna)
                                translated_seq = s_check2.translate()
                                if "*" in translated_seq:
                                    count3 = translated_seq.count("*")
                                    if min_num_index([count1, count2, count3]) == 0:
                                        pass
                                    elif min_num_index([count1, count2, count3]) == 1:
                                        addN = addN + 1
                                    elif min_num_index([count1, count2, count3]) == 2:
                                        addN = addN + 2
                                else:
                                    addN = addN + 2
                            else:
                                addN = addN + 1
                        else:
                            pass

                        if addN != 0:
                            extNs_new = (extNs + addN)%3
                            if extNs_new == 3:
                                extNs_new = 0
                            
                            record[i].seq = Seq("N"*extNs_new, generic_dna) + Seq(str(record[i].seq[extNs:]), generic_dna)
                            extNR = 3 - len(record[i].seq)%3
                            record[i].seq = record[i].seq + Seq(extNR*"N", generic_dna)

                        if record[i].seq[:3].count("N") == 3:
                            record[i].seq = record[i].seq[3:]
                
                        record[i].seq = ungap(record[i].seq)
                                
                    else:
                        record[i].seq = Seq(left_Ns, generic_dna) + Seq("".join([x for x in str(rec.seq) if x != "-"]), generic_dna) + Seq(right_Ns, generic_dna)
                        
                        if record[i].seq[:3].count("N") == 3:
                            record[i].seq = record[i].seq[3:]
                        
                        record[i].seq = ungap(record[i].seq)
                else:
                    record[i].seq = Seq(left_Ns, generic_dna) + Seq("".join([x for x in str(rec.seq) if x != "-"]), generic_dna) + Seq(right_Ns, generic_dna)
                    
                    if record[i].seq[:3].count("N") == 3:
                        record[i].seq = record[i].seq[3:]
                    
                    record[i].seq = ungap(record[i].seq)

                N_dict[rec.id] = (3 - record[i].seq[:3].count("N"))
                if record[i].seq[:3].count("N") == 3:
                    record[i].seq = record[i].seq[3:]

            SeqIO.write(record, fp, "fasta")

        store_post_Ns = dict()
        for rec in record:
            if "Homo_sapiens" in rec.id:
                pos_left_Ns = rec.seq[:3].count("N")

        for i, rec in enumerate(record):
            taxa_left_Ns = rec.seq[:3].count("N")
            if taxa_left_Ns != pos_left_Ns:
                store_post_Ns[rec.id] = taxa_left_Ns - pos_left_Ns

        cdsAlign(filename[:-4] + "withN.fas", filename[:-4] + "_aligned.fas", stage="mapper1")
        handle = open(filename[:-4] + "_aligned.fas", "r")
        record = list(SeqIO.parse(handle, "fasta"))
        for i, rec in enumerate(record):
            if "Homo_sapiens" in rec.id:

                if right_Ns == "":
                    pos_e = -1
                else:
                    pos_e = rec.seq[::-1].find(right_Ns) + len(right_Ns) - 1
                
                if left_Ns == "":
                    pos_s = -1
                else:
                    pos_s = rec.seq.find(left_Ns) + len(left_Ns) - 1
                
                break


        for i, rec in enumerate(record):

            if pos_e == -1:
                record[i].seq = rec.seq[pos_s+1:]
            else:
                record[i].seq = rec.seq[pos_s+1:-(pos_e+1)]

        with open(filename[:-4] + "_aligned.fas", "w") as fp:
            SeqIO.write(record, fp, "fasta")

        shutil.copy2(filename[:-4] + "_aligned.fas", geneName + "/" + filename[:-4].split("/")[-1] + ".fas")
        shutil.copy2(filename[:-4] + "withN.fas",  geneName + "/" + filename[:-4].split("/")[-1] + "_withN.fas")
        handle = open(filename[:-4] + "_aligned.fas", "rU")
        record = list(SeqIO.parse(handle, "fasta"))
        for rec in record:
            if "Homo_sapiens" in rec.id:
                n_left = rec.seq.strip("-")[:3].count("N")
                n_right = rec.seq.strip("-")[-3:].count("N")
                positions = [i for i, x in enumerate(rec.seq) if x == "-"]

        seqRecordObj = list()
        for rec in record:
            sequenceObj = Seq("", generic_dna)
            for i, nuc in enumerate(rec.seq):
                if i not in positions:
                    sequenceObj = sequenceObj + Seq(str(nuc), generic_dna)

            if n_right != 0:
                sequenceObj = sequenceObj[n_left:-n_right]
            else:
                sequenceObj = sequenceObj[n_left:]

            if len([x for x in sequenceObj if x != "N" and x != "-"]) != 0:
                if sequenceObj[:3].count("N") == 3:
                    sequenceObj = sequenceObj[3:]
                seqRecordObj.append(SeqRecord(Seq(str(sequenceObj), generic_dna), id=rec.id, name=rec.name, description=rec.description))

        for rec in seqRecordObj:
            if "Homo_sapiens" in rec.id:
                hum_len = len(rec.seq)

        for i, rec in enumerate(seqRecordObj):
            if len(rec.seq) > hum_len:
                seqRecordObj[i].seq = rec.seq[:hum_len]
            elif len(rec.seq) < hum_len:
                seqRecordObj[i].seq = rec.seq + Seq((hum_len - len(rec.seq))*"N", generic_dna)


        with open(filename[:-4] + ".nex", "w") as fp:
            SeqIO.write(seqRecordObj, fp, "nexus")
            os.remove(filename[:-4] + ".fas")
            os.remove(filename[:-4] + "_aligned.fas")
            os.remove(filename[:-4] + "withN.fas")
            try:
                pdist_values, mean_val, sd_val, pdist_values_nuc, mean_val_nuc, sd_val_nuc, cutoff_aa, cutoff_nuc = check_pdist_all(seqRecordObj, left_Ns, right_Ns)


                fp.write("\nbegin ConCat_Bin;\n\n")
            
                fp.write("\t[Amino acid sequence p distance stats]\n\n")
                store_outliers = list()
                for key, val in pdist_values.items():
                
                    if val < cutoff_aa:
                        if ("Monodelphis" in key or "Sarcophilus" in key) and val > 0.75:
                            continue
                        
                        message_writer = "-> Potential outlier"
                        store_outliers.append(key[0])
                        store_outliers.append(key[1])
                    else:
                        message_writer = ""
                
                    fp.write("\t%s:\t%s\t%s\n" %(key, val, message_writer))

                try:
                    frequency_dict = collect_frequency(store_outliers)
                    outliers = check_outlier(frequency_dict)
                except:
                    outliers = []


                if outliers != []:
                    aa_status = True
                else:
                    aa_status = False

                fp.write("\n\tMean:\t%s\n" %mean_val)
                fp.write("\tStandard Deviation:\t%s\n" %sd_val)
                fp.write("\tOutliers: \t%s\n\n\n" %outliers)
                fp.write("\n\n\t[Nucleotide sequence p distance stats]\n\n")
                for key, val in pdist_values_nuc.items():
                    if val < cutoff_nuc:
                        if ("Monodelphis" in key or "Sarcophilus" in key) and val > 0.75:
                                continue
                    
                        message_writer_nuc = "-> Potential outlier"
                        store_outliers.append(key[0])
                        store_outliers.append(key[1])
                    else:
                        message_writer_nuc = ""

                    fp.write("\t%s:\t%s\t%s\n" %(key, val, message_writer_nuc))

                try:
                    frequency_dict = collect_frequency(store_outliers)
                    outliers_nuc = check_outlier(frequency_dict)
                except:
                    outliers_nuc = None

                if outliers_nuc != []:
                    nuc_status = True
                else:
                    nuc_status = False

                if aa_status == True:
                    for taxaname in outliers:
                        outlier_files.write("%s\t%s\t%s\n" %(geneName, filename.split("/")[-1], taxaname))

                fp.write("\n\tMean:\t%s\n" %mean_val_nuc)
                fp.write("\tStandard Deviation:\t%s\n" %sd_val_nuc)
                fp.write("\tOutliers: \t%s\n\n" %outliers)
                fp.write("end;\n")
                    
            except:
                pass
                


def collect_frequency(listObj):
    retDict = dict()
    for x in listObj:
        if x not in retDict.keys():
            retDict[x] = 1
        else:
            retDict[x] = retDict[x] + 1

    return retDict


def check_outlier(dictObj):
    retList = list()
    fval = [val for key, val in dictObj.items()]
    fmean = mean(fval)
    fstdev = stddev(fval)
    for key, val in dictObj.items():
        if val > fmean + 2*fstdev:
            retList.append(key)

    return retList


def check_pdist_all(recordObj_final, left_Ns, right_Ns):
    potential_outliers_amino = list()
    potential_outliers_nuc = list()
    seq_ids = [x.id for x in recordObj_final]
    store_used_id = list()
    pdist_values = dict()
    pdist_values_nuc = dict()
    for rec in recordObj_final:
        store_used_id.append(rec.id)
        for rec_in in recordObj_final:
            if rec_in.id not in store_used_id:
                seq1 = Seq(left_Ns, generic_dna) + rec.seq + Seq(right_Ns, generic_dna)
                seq2 = Seq(left_Ns, generic_dna) + rec_in.seq + Seq(right_Ns, generic_dna)
                pdist_values[(rec.id, rec_in.id)] = pdist(Seq(str(seq1).replace("-", "N"), generic_dna).translate(), Seq(str(seq2).replace("-", "N"), generic_dna).translate())
                pdist_values_nuc[(rec.id, rec_in.id)] = pdist(rec.seq, rec_in.seq)

    pdist_vals = [val for key, val in pdist_values.items()]
    pdist_vals_nuc = [val for key, val in pdist_values_nuc.items()]
    mean_val = mean(pdist_vals)
    mean_val_nuc = mean(pdist_vals_nuc)
    sd_val = stddev(pdist_vals)
    sd_val_nuc = stddev(pdist_vals_nuc)
    cutoff_aa = mean_val - (2*sd_val)
    cutoff_nuc = mean_val - (2*sd_val_nuc)
    if cutoff_aa > 0.9:
        cutoff_aa = 0.9
    if cutoff_nuc > 0.9:
        cutoff_nuc = 0.9
                                                
    return pdist_values, mean_val, sd_val, pdist_values_nuc, mean_val_nuc, sd_val_nuc, cutoff_aa, cutoff_nuc


            
def concatenate_exons(listObject_human, output1):
    for fileobj in listObject_human:
        file_list = glob.glob(output1+"/"+fileobj+"/*.nex")
        nexi = [(fname, Nexus.Nexus(fname)) for fname in natural_sort(file_list)]
        combined = Nexus.combine(nexi)
        with open(output1+"/"+fileobj+".nex", 'w') as fp:
            combined.write_nexus_data(fp)


####################################################################################
# BEGIN EXECUTION
####################################################################################


def execute_map(listObject_h, output1):
    try:
        os.mkdir(argmnts.o)
    except:
        pass

    logObj = "exonName.log"

    nContentRet, listObject_human, data_dict_human, record_original_tot = exec_mapping(listObject_h, tag="human")
    exon_failed = align_exon(nContentRet, listObject_h)

    if exon_failed:
        print("Failed to import sequences for following exons:\n%s" %exon_failed)

    exon_file_names = dict()

    for geneName in listObject_h:
        try:
            exon_file_names[geneName] = exon_names(logObj, data_dict_human, record_original_tot, geneName)
        except:
            pass


    outlier_files = open("Outliers.log", "w")
    outlier_files.write("Gene\tExon\tAmino Acid Outlier\n\n")

    for geneName in listObject_h:
        try:
            os.mkdir(geneName)
        except:
            pass

        try:
            mapper_whole(geneName, nContentRet, exon_file_names[geneName], outlier_files)
            outlier_files.write("\n\n")
        except KeyError:
            continue

    outlier_files.close()

    transfer_to(output1, listObject_h)
    concatenate_exons(listObject_human, output1)




parser = argparse.ArgumentParser(prog='Consensus',
                                 version= 'Consensus-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t\t Welcome to Consensus-1.0 alignment handling program
    \t\t\t ConCat-Align conducts sequence alignment using Muscle/Mafft
    \t\t\t sequence alignment programs.
    \t\t\tWritten by Ambuj Kumar, University of Florida
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))


parser.add_argument('-i', type=str, default=None,
                    help='gene names')

parser.add_argument('-o', type=str, default=None,
                    help='CDS alignemnt output folder name')


argmnts = parser.parse_args()

if argmnts.i == None:
    parser.error('-i argument is required.')

if argmnts.o == None:
    argmnts.o = "Output"

try:
    os.mkdir(argmnts.o)
except OSError:
    shutil.rmtree(argmnts.o)
    os.mkdir(argmnts.o)


def main():
    try:
        listObject_h = [x for x in open(argmnts.i, 'r').readlines() if x != '' and x != '\n']
    except IOError:
        print("\n%s file not found." %argmnts.i)

    execute_map(listObject_h, argmnts.o)


if __name__ == "__main__":
    main()




