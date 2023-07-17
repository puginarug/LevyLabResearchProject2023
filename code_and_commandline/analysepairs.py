#!/usr/bin/env python3

import sys
import re
import numpy as np
import argparse
from Bio.PDB import PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser

def pDockQ(structure,cutoff,if1,if2,disocut1=50,disocut2=70,disocut3=90):
    n=0
    tiny=1.e-20
    atoms=[]
    coords=[]
    residues=[]
    chains=[]
    for chain in structure.get_chains():
        chains+=[chain]
        n+=1

        #print (chains)
    interface_residues1=[]
    interface_residues2=[]
    if n>1:
        for res1 in chains[0]:
            for res2 in chains[1]:
                #Atom-atom distance
                #print (res1,res2)
                test=False
                for i in res1:
                    if test:break
                    for j in res2:
                        dist = np.linalg.norm(i.coord - j.coord)
                        #scipy.spatial.distance.euclidian
                        #dist = distance.euclidean(coords[i], coords[len(ch1_res_nos)+j]) #Need to add l1 to get the right coords
                        if dist < cutoff:
                            #Save residues
                            #print ("Appending",res1,res2)
                            interface_residues1.append(res1.id[1])
                            interface_residues2.append(res2.id[1])
                            test=True
                            break
                        elif dist > 2*cutoff: # To speed up things
                            test=True
                            break

    #print (interface_residues1,interface_residues2)
    #print (np.unique(interface_residues1),np.unique(interface_residues2))
    #interface_res_num.append(np.unique(interface_residues).shape[0])
    #atoms, residue_numbers, coords = np.array(atoms), np.array(residue_numbers), np.array(coords)
    #print("test2",if1.shape[0],if2.shape[0])
    if (if1.size==1):
        if1=np.unique(interface_residues1)
    if (if2.size==1):
        if2=np.unique(interface_residues2)
    NumRes=if1.shape[0]+if2.shape[0]
    i=tiny
    b=0
    b1=0
    b2=0
    i1=0
    i2=0
    NumDiso1=[0,0,0,0]
    NumDiso2=[0,0,0,0]
    NumIFDiso1=0
    NumIFDiso2=0
    for res in chains[0]:
        b1+=res['CA'].get_bfactor()
        i1+=1
        if res['CA'].get_bfactor()>disocut3: # >90
            NumDiso1[0]+=1
        elif res['CA'].get_bfactor()>disocut2: # 70-90
            NumDiso1[1]+=1
        elif res['CA'].get_bfactor()>disocut1: # 50-70
            NumDiso1[2]+=1
        else: # <50
            NumDiso1[3]+=1
        if res.id[1] in if1:
            b+=res['CA'].get_bfactor()
            i+=1
            if res['CA'].get_bfactor()<disocut1: # >90
                NumIFDiso1+=1
    if n>1:
        for res in chains[1]:
            b2+=res['CA'].get_bfactor()
            i2+=1
            if res['CA'].get_bfactor()>disocut3: # >90
                NumDiso2[0]+=1
            elif res['CA'].get_bfactor()>disocut2: # 70-90
                NumDiso2[1]+=1
            elif res['CA'].get_bfactor()>disocut1: # 50-70
                NumDiso2[2]+=1
            else: # <50
                NumDiso2[3]+=1
            if res.id[1] in if2:
                b+=res['CA'].get_bfactor()
                i+=1
                if res['CA'].get_bfactor()<disocut1: # >90
                    NumIFDiso2+=1
    else:
        b2=b1
        i2=i1
        NumDiso2=NumDiso1
    IF_plDDT=b/i
    plDDT1=b1/i1
    plDDT2=b2/i2
    #print ("test",b,b1,b2,i,i1,i2,NumDiso1,NumDiso2)
    #Get res nos
    #Get chain cut
    #ch1_res_nos = np.argwhere(residue_numbers<=l1)[:,0] #All residue numbers
    #ch2_res_nos =  np.argwhere(residue_numbers>l1)[:,0]
    
    #print (NumRes,IF_plDDT)
    return (NumRes,IF_plDDT,plDDT1,plDDT2,NumDiso1,NumDiso2,i1,i2,if1,if2,NumIFDiso1,NumIFDiso2)



def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return (y)

popt=[7.07140240e-01, 3.88062162e+02, 3.14767156e-02, 3.13182907e-02]

arg_parser = argparse.ArgumentParser(description="Calculates pDockQ from NumRes and IF_plDDT")
group = arg_parser.add_mutually_exclusive_group(required=True)
group.add_argument("-p","--pdb",type=argparse.FileType('r'),help="Input pdb file of complex with pLddt values in bfactor columns")
group.add_argument("-c","--cif",type=argparse.FileType('r'),help="Input cif file of complex with plddt values in bfactor columns")
group1 = arg_parser.add_mutually_exclusive_group(required=True)
group1.add_argument("-p1","--pdb1",type=argparse.FileType('r'),help="Input pdb file of chain 1 with pLddt values in bfactor columns")
group1.add_argument("-c1","--cif1",type=argparse.FileType('r'),help="Input cif file of chain 1 with plddt values in bfactor columns")
group2 = arg_parser.add_mutually_exclusive_group(required=True)
group2.add_argument("-p2","--pdb2",type=argparse.FileType('r'),help="Input pdb of chain2 with file pLddt values in bfactor columns")
group2.add_argument("-c2","--cif2",type=argparse.FileType('r'),help="Input cif of chain 2 file plddt values in bfactor columns")

arg_parser.add_argument("-v","--verbose", action='store_true', required=False,help="Average plDDT in interface residues")

args = arg_parser.parse_args()
cutoff=10

if (args.cif):
    Name=args.cif.name
    bio_parser = MMCIFParser()
    structure_file = args.cif
    structure_id = args.cif.name[:-4]
    structure = bio_parser.get_structure(structure_id, structure_file)
elif (args.pdb):
    Name=args.pdb.name
    bio_parser = PDBParser()
    structure_file = args.pdb
    structure_id = args.pdb.name[:-4]
    structure = bio_parser.get_structure(structure_id, structure_file)
    
if (args.cif1):
    bio_parser = MMCIFParser()
    structure1_file = args.cif1
    structure1_id = args.cif1.name[:-4]
    structure1 = bio_parser.get_structure(structure1_id, structure1_file)
elif (args.pdb1):
    bio_parser = PDBParser()
    structure1_file = args.pdb1
    structure1_id = args.pdb1.name[:-4]
    structure1 = bio_parser.get_structure(structure1_id, structure1_file)
    
if (args.cif2):
    bio_parser = MMCIFParser()
    structure2_file = args.cif2
    structure2_id = args.cif2.name[:-4]
    structure2 = bio_parser.get_structure(structure2_id, structure2_file)
elif (args.pdb2):
    bio_parser = PDBParser()
    structure2_file = args.pdb2
    structure2_id = args.pdb2.name[:-4]
    structure2 = bio_parser.get_structure(structure2_id, structure2_file)
    
tiny=1.e-20
cutoff2=3
#print (NumRes,tiny,IF_plDDT, popt)


if1=np.empty([])
if2=np.empty([])
NumRes,IF_plDDT,plDDT1,plDDT2,Diso1,Diso2,len1,len2,if1,if2,NumIFDiso1,NumIFDiso2=pDockQ(structure,cutoff,if1,if2)
NumResOverlap,IF_plDDTOverlap,plDDT1Overlap,plDDT2overlap,Diso1overlap,Diso2overlap,len1,len2,if1overlap,if2overlap,NumIFDiso1overlap,NumIFDiso2overlap=pDockQ(structure,cutoff2,np.empty([]),np.empty([]))
        #print (NumRes,IF_plDDT,plDDT1,plDDT2,Diso1,Diso2,len1,len2)
pd=sigmoid(np.log(NumRes+tiny)*IF_plDDT,*popt)
#print (pd)
NumRes1,IF_plDDT1,plDDT11,plDDT21,Diso11,Diso21,len11,len21,if1,temp,NumIFDiso11,temp2=pDockQ(structure1,cutoff,if1,np.empty([]))
#print (NumRes,IF_plDDT,plDDT1,plDDT2,Diso1,Diso2,len1,len2)

NumRes2,IF_plDDT2,plDDT12,plDDT22,Diso12,Diso22,len12,len22,if2,temp,NumIFDiso12,temp2=pDockQ(structure2,cutoff,if2,np.empty([]))
#print (NumRes,IF_plDDT,plDDT1,plDDT2,Diso1,Diso2,len1,len2)

Name=re.sub(r'.*/','',Name)
Name=re.sub(r'.pdb$','',Name)
Name=re.sub(r'.cif$','',Name)

print ("Name,pDockQ,NumRes,IF_plDDT,IF_plDDT1,IF_plDDT2,plDDT1,plDDT2,NumDiso1,NumDiso2,NumOverlap,len1,len2,plDDT1s,plDDT2s,NumDiso1s,NumDiso2s,NumIFDiso1,NumIFDiso2,NumIFDiso1monomer,NumIFDiso2monomer" )

print ("%s,%f,%d,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%f,%f,%d,%d,%d,%d,%d,%d" %( Name,pd,NumRes,IF_plDDT,IF_plDDT1,IF_plDDT2,plDDT1,plDDT2,Diso1[3],Diso2[3],NumResOverlap,len1,len2,plDDT11,plDDT12,Diso11[3],Diso12[3],NumIFDiso1,NumIFDiso2,NumIFDiso11,NumIFDiso12 ))

    
# NumRes...
