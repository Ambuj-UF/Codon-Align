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
import glob
import shutil
import platform
import subprocess

from src.aligner import cdsAlign

from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Align import MultipleSeqAlignment






def to_list(str_obj): return [x for x in str_obj.lstrip("[").rstrip("]").split(", ")]

def reverser(a1, a2): return (a1, a2)[1], (a1, a2)[0]

def file_is_empty(path): return os.stat(path).st_size==0

def pick_longest(inputObject): return [rec for rec in inputObject if len(rec.seq) >= max([len(x.seq) for x in inputObject])][0]

def strfind(s, ch): return [i for i, ltr in enumerate(s) if ltr == ch]

def pdist(seq1, seq2): return float(sum([1 for n1, n2 in zip(seq1, seq2) if str(n1) == str(n2)]))/len(seq1)

def key_max(dictObj): return max(dictObj.iterkeys(), key=lambda k: dictObj[k])

def combine_dict(a, b): return dict(a.items() + b.items() + [(k, a[k] + b[k]) for k in set(b) & set(a)])

def pdist(s1, s2): return float(len([x for x, y in zip(s1, s2) if x == y and x != "-" and y != "-"]))/len([x for x, y in zip(s1, s2) if x != "-" and y != "-"])


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def transfer_to(output):
    try:
        os.mkdir(output)
    except:
        pass

    folders = glob.glob("Sequences/*/Aligned/")
    #folders = glob.glob("Sequences/ASPM/Aligned/")
    for folder in folders:
        recursive_overwrite(folder, output + "/" + folder.split("/")[1])

    print("Files are in ConCat data directory\n")



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



def exon_names(logfileName):
    file_folders = glob.glob("Sequences/*/Aligned/")
    #file_folders = glob.glob("Sequences/MCPH1/Aligned/")
    #file_folders = glob.glob("Sequences/ASPM/Aligned/")
    logData = open(logfileName, 'a+')

    for folders in file_folders:
        print folders
        files = natural_sort(glob.glob(folders + "*.*"))
        for filename in files:
            handle = open(filename, "rU")
            record = list(SeqIO.parse(handle, "nexus"))
            handle.close()
            for i, rec in enumerate(record):
                if "gi|" in rec.id:
                    record[i].id = "Homo_sapiens_CCDS|CCDS"
        
            with open(filename, "w") as fp:
                SeqIO.write(record, fp, "nexus")
    
        handle = open(files[0], 'rU')
        record = list(SeqIO.parse(handle, "nexus"))
        handle.close()
        for rec in record:
            if "Homo_sapiens_CCDS" in rec.id and str(rec.seq[:3]) == "ATG":
                files = files
                break
            elif "Homo_sapiens_CCDS" in rec.id and str(rec.seq[:3]) != "ATG":
                files = files[::-1]
                break
    
        for i, filename in enumerate(files):
            os.rename(filename, "/".join(filename.split("/")[:-1]) + "/exon" + str(i+1) + ".nex")
            logData.write("%s - exon%s"%(filename, i))

    logData.close()


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

            while flagThing != True:
                try:
                    tblastn_cline = NcbitblastnCommandline(cmd='tblastn',
                                              db="%s"%(blast_to),
                                              num_threads=6,
                                              outfmt=5,
                                              seg="no",
                                              out=("%s/Result%s.xml"%(strore_blast_res,i))
                                              )
                    
                    flagThing = True
                except (urllib2.URLError, httplib.HTTPException), err:
                    continue


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
                        #fcount = fcount + 1
                        #if fcount > 1:
                        #    break
                                
                    if hitFlag == True:
                        if re.search('<Hit_def>', line) != None:
                            description = line.strip().strip("\n")[9:-10]
                            id = line.strip().strip("\n")[9:-10]
                            retData[rec.id][id] = dict()
                            #description = line.strip().rstrip().strip("<Hit_def>").strip('</')
                            retData[rec.id][id]["description"] = description
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
                    handle_in = Entrez.efetch(db="nucleotide",
                                       id=ccds_object["id"],
                                       rettype="fasta",
                                       strand=+1,
                                       seq_start=seq_coord[0] + 1,
                                       seq_stop=seq_coord[1] + 1)
                                       
                    record = SeqIO.read(handle_in, "fasta")
                    flagThing = True
        
                except (IOError, httplib.HTTPException):
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

        
            print (geneName, record.id, record.seq)
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
                #print inkey, inval["eval"]
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
        
            #if x.seq[-1:] == "X":
            #    if i == len(ccdsObj) - 1:
            #        newSeq = ccdsObj[i-1].seq + x.seq
            #        idObj = x.id + "+" + ccdsObj[i-1].id + "_lag"
            #    elif i == 0:
            #        newSeq = ccdsObj[i+1].seq + x.seq
            #        idObj = x.id + "+" + ccdsObj[i+1].id + "_lag"
            #    else:
            #        newSeq = ccdsObj[i+1].seq + x.seq
            #        idObj = x.id + "+" + ccdsObj[i+1].id + "_lag"
            #
            #elif x.seq[:1] == "M":
            #    if i == len(ccdsObj) - 1:
            #        newSeq = x.seq + ccdsObj[i-1].seq
            #        idObj = x.id + "+" + ccdsObj[i-1].id + "_lead"
            #    elif i == 0:
            #        newSeq = x.seq + ccdsObj[i+1].seq
            #        idObj = x.id + "+" + ccdsObj[i+1].id + "_lead"
            #    else:
            #        newSeq = x.seq + ccdsObj[i+1].seq
            #        idObj = x.id + "+" + ccdsObj[i+1].id + "_lead"
            #else:
            #    if i == len(ccdsObj) - 1:
            #        newSeq = ccdsObj[i-1].seq + x.seq
            #        idObj = x.id + "+" + ccdsObj[i-1].id + "_lag"
            #    else:
            #        newSeq = x.seq + ccdsObj[i+1].seq
            #        idObj = x.id + "+" + ccdsObj[i+1].id + "_lead"

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
                    print leadElement
                    print prefiles
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
                                #recordObj_twins[i].seq = rec.seq[:-len(leadRecord[rec.id.split("|")[0]])]
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
                                #recordObj_twins[i].seq = rec.seq[len(lagRecord[rec.id.split("|")[0]]):]
                                s1 = str(lagRecord[rec.id.split("|")[0]].seq[-10:]).strip("N")
                                s = str(rec.seq)
                                try:
                                    pos = re.search(str(s1), str(s)).end()
                                    recordObj_twins[i].seq = rec.seq[pos:]
                                    #recordObj_twins[i].seq = rec.seq[re.search(str(lagRecord[rec.id.split("|")[0]][-10:]).strip("N"), str(rec.seq)).end():]
                                except AttributeError:
                                    recordObj_twins[i].seq = rec.seq[len(lagRecord[rec.id.split("|")[0]]):]
                                    continue
                            except KeyError:
                                continue


                newRecObject_twins = list()
                for rec in recordObj_twins:
                    if len(rec.seq) != 0:
                        newRecObject_twins.append(rec)

                with open(filename.split("+")[0].split("gi")[0] + filename.split("+")[0].split(":")[-1] + ".fas", "w") as fp:
                    os.remove(filename)
                    SeqIO.write(recordObj_twins, fp, "fasta")
                    
                        
            
    return nContentRet, listObject_human, data_dict_human, record_original_tot


##########################################################################################################
# Exon alignments
##########################################################################################################


def align_exon(nContentRet):
    exon_failed = list()
    seqFolders = glob.glob("Sequences/*/")
    #seqFolders = ["Sequences/ASPM/"]

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

            with open(filename, "w") as fp:
                for i, rec in enumerate(record):
                    record[i].seq = Seq(left_Ns, generic_dna) + rec.seq + Seq(right_Ns, generic_dna)
            
                SeqIO.write(record, fp, "fasta")

            outfile = filename.replace("Results", "Aligned")
            print filename
            
            try:
                cdsAlign(filename, outfile)
            except:
                exon_failed.append(filename)
                with open(outfile[:-3] + "nex", "w") as fp:
                    SeqIO.write(record, fp, "nexus")
            
                continue
            
            handle = open(outfile, 'rU')
            record = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            for rec in record:
                if "gi|" in rec.id:
                    positions = [i for i, x in enumerate(rec.seq) if x == "-"]
    
            seqRecordObj = list()
            for rec in record:
                sequenceObj = Seq("", generic_dna)
                for i, nuc in enumerate(rec.seq):
                    if i not in positions:
                        sequenceObj = sequenceObj + Seq(str(nuc), generic_dna)
            
                if len([x for x in rec.seq if x != "N" and x != "-"]) != 0:
                    seqRecordObj.append(SeqRecord(Seq(str(sequenceObj), generic_dna), id=rec.id, name=rec.name, description=rec.description))
                    #print Seq(str(rec.seq), generic_dna)
                    
            ref_seq_record = list()
            human_seq = [x.seq for x in seqRecordObj if "gi|" in x.id][0]
            with open("pdist_log.txt", "a+") as fp:
                fp.write("\n%s\n\n"%folders.split("/")[-1])
                for rec in seqRecordObj:
                    p_distance = pdist(rec.seq, human_seq)
                    if p_distance >= 0.5:
                        ref_seq_record.append(rec)
                    
                    fp.write("%s\t%s\t%s\n" %(filename, rec.id, p_distance))
            

        
            with open(outfile[:-3] + "nex", "w") as fp:
                SeqIO.write(ref_seq_record, fp, "nexus")

            os.remove(outfile)

    logfile.close()

    return exon_failed

####################################################################################
# BEGIN EXECUTION
####################################################################################

#genes = open("cong_micro.txt", 'r').readlines()
#listObject_h = [x.rstrip("\n") for x in genes]


def fetcher(listObject_h, output, logObj = "exonName.log"):
    nContentRet, listObject_human, data_dict_human, record_original_tot = exec_mapping(listObject_h, tag="human")
    exon_failed = align_exon(nContentRet)
    print("Failed to import sequences for following exons:\n%s" %exon_failed)
    exon_names(logObj)
    transfer_to(output)








