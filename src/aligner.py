################################################################################################################
# Alignment editing program. It removes indels and stop codons and provides an option to adjust the alignment  #
# according to user defined protein or species sequence                                                        #
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


import sys
import os
import platform
import subprocess


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC
from Bio.Seq import translate
from Bio.Alphabet import SingleLetterAlphabet



def manage_seqLength(x): return max(zip((x.count(item) for item in set(x)), set(x)))[1]

def _spliter(str, num): return [ str[start:start+num] for start in range(0, len(str), num) ]


def _groupy(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def _translator(recordData, ign, omit, table):
    proteinSeqList = list()
    recordsFunc = recordData
    
    for i, rec in enumerate(recordsFunc):
        counter = dict()
        if len(rec.seq)%3 != 0:
            rec.seq = rec.seq[:-(len(rec.seq)%3)]

        seqT = rec.seq.translate()
        
        for j, obj in enumerate(seqT):
            if '*' in obj:
                seqT = seqT[:j] + 'Z' + seqT[j+1:]
        
        proteinSeqList.append(SeqRecord(Seq(str(seqT), IUPAC.protein), id=rec.id, name=rec.name, description=rec.description))


    with open('translated.fas', 'w') as fp:
        SeqIO.write(proteinSeqList, fp, 'fasta')

    
    return recordsFunc




def _alignP(pkg, arguments=None):
    if pkg == 'muscle':
        if 'Darwin' in platform.system():
            subprocess.call("./src/muscle/muscle -in translated.fas -out tAligned.fas", shell=True)
        else:
            subprocess.call("./src/muscle/muscleLinux -in translated.fas -out tAligned.fas", shell=True)
    
    else:
        arguments = arguments.replace('[', '').replace(']', '')
        subprocess.call("./src/mafft/mafft.bat %s translated.fas > tAligned.fas" %arguments, shell=True)


def _cleanAli(recordNuc, omit, fileName, stage):
    handleP = open('tAligned.fas', 'rU')
    records = list(SeqIO.parse(handleP, 'fasta'))
    
    store = list()
    two_state_flag = False
    
    for rec in recordNuc:
        if "gi|" in rec.id or "Homo_sapiens" in rec.id:
            n_count_s = rec.seq[:3].count("N")
            n_count_e = rec.seq[-3:].count("N")
            break
        elif "Querry" in rec.id:
            n_count_s = 0
            n_count_e = 0
            two_state_flag = True


    for i, rec in enumerate(records):
        
        nucData = [x.seq for x in recordNuc if x.id in rec.id]
        nucSeqData = _spliter(nucData[0], 3)
        sequence = Seq("", SingleLetterAlphabet()); pos = 0
        
        #print len([x for x in rec.seq if x!="-"]), len(nucSeqData)
        for j, amino in enumerate(rec.seq):
            if amino == '-':
                sequence = sequence + Seq("---", SingleLetterAlphabet())
            elif amino == "Z":
                sequence = sequence + Seq("NNN", SingleLetterAlphabet())
                pos = pos + 1
            else:
                if pos == 0 or pos == len(nucSeqData) - 1:
                    if pos == 0:
                        tot_N_s = nucSeqData[pos].count("N")
                    
                    if stage == "mapper":
                        sequence = sequence + nucSeqData[pos]
                    else:
                        sequence = sequence + nucSeqData[pos].strip("N")
                else:
                    sequence = sequence + nucSeqData[pos]
                
                pos = pos + 1


        if two_state_flag == False:
            if stage == "mapper":
                records[i].seq = Seq(str(sequence), SingleLetterAlphabet())
                
                #if n_count_e != 0:
                #    records[i].seq = Seq("N"*(tot_N_s-n_count_s) + str(sequence)[n_count_s:-n_count_e], SingleLetterAlphabet())
                #else:
                #    records[i].seq = Seq("N"*(tot_N_s-n_count_s) + str(sequence)[n_count_s:], SingleLetterAlphabet())
            else:
                records[i].seq = Seq("N"*(tot_N_s-n_count_s) + str(sequence).strip("N"), SingleLetterAlphabet())
        else:
            records[i].seq = Seq(str(sequence).strip("N"), SingleLetterAlphabet())
        
    #if stage != "mapper":
    #    optimal_length = manage_seqLength([len(rec.seq) for rec in records])

    #    for i, rec in enumerate(records):
    #        rec.seq = rec.seq[:optimal_length]

    with open(fileName, 'w') as fp:
        SeqIO.write(records, fp, "fasta")


    os.remove('translated.fas')
    os.remove('tAligned.fas')





def _cleanAli2(recordNuc, omit, fileName, stage):
    handleP = open('tAligned.fas', 'rU')
    records = list(SeqIO.parse(handleP, 'fasta'))
    
    store = list()
    
    for rec in records:
        if "gi|" in rec.id or "Homo_sapiens" in rec.id:
            n_count_s = rec.seq[:3].count("N")
            n_count_e = rec.seq[-3:].count("N")
            break

    #print records
    #print recordNuc

    for i, rec in enumerate(records):
        
        nucData = [x.seq for x in recordNuc if x.id in rec.id]
        nucSeqData = _spliter(nucData[0], 3)

        if stage == "mapper":
            nucSeqData[0] = nucSeqData[0].lstrip("N")
            nucSeqData[-1] = nucSeqData[-1].rstrip("N")
        sequence = Seq("", SingleLetterAlphabet())
        pos = 0
        
        for j, amino in enumerate(rec.seq):
            if amino == '-':
                sequence = sequence + Seq("---", SingleLetterAlphabet())
            elif amino == "Z":
                sequence = sequence + Seq("NNN", SingleLetterAlphabet())
                pos = pos + 1
            else:
                sequence = sequence + nucSeqData[pos]
                pos = pos + 1
        
        records[i].seq = Seq(str(sequence), SingleLetterAlphabet())



    with open(fileName, 'w') as fp:
        SeqIO.write(records, fp, "fasta")
    
    
    os.remove('translated.fas')
    os.remove('tAligned.fas')







def cdsAlign(inputFile, outFile, pkg='muscle', ign=True, CT=None, stage="noMapping"):
    
    
    
    codonTables = ['Ascidian Mitochondrial', 'SGC9', 'Coelenterate Mitochondrial', 'Protozoan Mitochondrial', 'Vertebrate Mitochondrial', 'Plant Plastid', 'Thraustochytrium Mitochondrial', 'Blepharisma Macronuclear', 'Mold Mitochondrial', 'Invertebrate Mitochondrial', 'Standard', 'Trematode Mitochondrial', 'Scenedesmus obliquus Mitochondrial', 'Euplotid Nuclear', 'Yeast Mitochondrial', 'Spiroplasma', 'Alternative Flatworm Mitochondrial', 'Ciliate Nuclear', 'SGC8', 'Alternative Yeast Nuclear', 'Hexamita Nuclear', 'SGC5', 'SGC4', 'SGC3', 'SGC2', 'SGC1', 'SGC0', 'Flatworm Mitochondrial', 'Dasycladacean Nuclear', 'Chlorophycean Mitochondrial', 'Mycoplasma', 'Bacterial', 'Echinoderm Mitochondrial']
    
    omit = False
    
    if CT == None:
        table = CodonTable.ambiguous_dna_by_id[1]
    elif CT != None and CT in codonTables:
        table = CodonTable.ambiguous_generic_by_name[CT]
    else:
        table = CodonTable.ambiguous_generic_by_name['Standard']

    
    handle = open(inputFile, 'rU')
    records = list(SeqIO.parse(handle, "fasta"))

    for j, rec in enumerate(records):
        #print rec.seq
        if 'TAA' in rec.seq[-3:] or 'TGA' in rec.seq[-3:] or 'TAG' in rec.seq[-3:]:
            records[j].seq = rec.seq[0:-3]



    records = _translator(records, ign, omit, table)
    _alignP(pkg)
    if stage == "mapper" or stage ==stage == "mapper1":
        _cleanAli2(records, omit, outFile, stage)
    else:
        _cleanAli(records, omit, outFile, stage)











