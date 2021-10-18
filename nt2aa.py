# Nucleotide to Amino Acid Translator
# Input file must be cluster count csv file with the following fields: gene,pos,base,A,C,G,T,N,
# These input files contain the total read counts for each nucleotide type
# Usage: python nt2aa.py -i inputfile.csv -o outputfile.csv -s start_pos_value -t threshold_aa_proportion_display
#   if the -v TRUE option is used the codon preference details are displayed in the output file  

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--treshold", help="file input name")
parser.add_argument("-s", "--start_position", help="date to filter data and get most cookie active")
parser.add_argument("-i", "--file_to_process", help="date to filter data and get most cookie active")
parser.add_argument("-o", "--output_file_name", help="date to filter data and get most cookie active")
parser.add_argument("-v", "--verbose", help="date to filter data and get most cookie active")

cases={"1":"2","2":"0","0":"1"}
args = parser.parse_args()
r=int(cases[str(int(args.start_position)%3)])
print(int(args.start_position)%3) 

cols="""gene
pos
sorted
AAA
AAC
AAG
AAT
ACA
ACC
ACG
ACT
AGA
AGC
AGG
AGT
ATA
ATC
ATG
ATT
CAA
CAC
CAG
CAT
CCA
CCC
CCG
CCT
CGA
CGC
CGG
CGT
CTA
CTC
CTG
CTT
GAA
GAC
GAG
GAT
GCA
GCC
GCG
GCT
GGA
GGC
GGG
GGT
GTA
GTC
GTG
GTT
TAA
TAC
TAG
TAT
TCA
TCC
TCG
TCT
TGA
TGC
TGG
TGT
TTA
TTC
TTG
TTT
A
C
D
E
F
G
H
I
K
L
M
N
P
Q
R
S
T
V
W
X
Y
ref
alt""".split("\n")
verbos="""gene
pos
AAA
AAC
AAG
AAT
ACA
ACC
ACG
ACT
AGA
AGC
AGG
AGT
ATA
ATC
ATG
ATT
CAA
CAC
CAG
CAT
CCA
CCC
CCG
CCT
CGA
CGC
CGG
CGT
CTA
CTC
CTG
CTT
GAA
GAC
GAG
GAT
GCA
GCC
GCG
GCT
GGA
GGC
GGG
GGT
GTA
GTC
GTG
GTT
TAA
TAC
TAG
TAT
TCA
TCC
TCG
TCT
TGA
TGC
TGG
TGT
TTA
TTC
TTG
TTT
A
C
D
E
F
G
H
I
K
L
M
N
P
Q
R
S
T
V
W
X
Y
ref
alt""".split("\n")
table = {
        "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
        "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
        "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
        "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
        "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
        "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
        "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
        "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
        "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
        "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
        "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
        "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
        "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
        "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
        "TAC":"Y", "TAT":"Y", "TAA":"X", "TAG":"X",
        "TGC":"C", "TGT":"C", "TGA":"X", "TGG":"W",
    }

import pandas
import os
mylist=" AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT, CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT, CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT, GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT, GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT, TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT, TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT"
mylist1="""A	= GCA + GCC + GCG + GCT
C	= TGC + TGT
D	= GAC + GAT
E	= GAA + GAG
F 	= TTC + TTT
G	= GGA + GGC + GGG + GGT
H 	= CAC + CAT
I 	= ATA + ATC + ATT
K 	= AAA + AAG
L 	= CTA + CTC + CTG + CTT + TTA + TTG
M 	= ATG
N 	= AAC + AAT
P 	= CCA + CCC + CCG + CCT
Q	= CAA + CAG
R 	= CGA + CGC + CGG + CGT + AGA + AGG
S 	= AGC + AGT + TCA + TCC + TCG + TCT
T 	= ACA + ACC + ACG + ACT
V 	= GTA + GTC + GTG + GTT
W 	= TGG
X	= TAA + TAG+ TGA
Y	= TAC + TAT"""
for u in [args.file_to_process]: #os.listdir("datas")
#topics_data = pd.DataFrame(new)

 p=0
 df = pandas.read_csv(u)
 #df=df[16:len(df)-int(len(df)/3)]
 df=df[df["pos"]>int(args.start_position)]
 
 pi=[]
 for pa in "ACGTN":
  df["P("+pa+")"]=(df[pa]/(df["A"]+df["C"]+df["G"]+df["T"]+df["N"])).round(3)
  pi.append("P("+pa+")")
  pi.append(pa)
 for comb in mylist.split(","):
  pi.append(comb[1:])   
  comb=comb[1:]
  print(comb)
 #df["P("+comb[0]+")"]=list(range(len(df["P("+comb[0]+")"])))
    #if not el=="A": 
    # df["P("+el+")"]=list(range(len(somme)))
  new=[]
  new.extend(df["P("+comb[1]+")"])
  new.extend([0]) 
  new1=[]
  new1.extend(df["P("+comb[2]+")"])    
  new1.extend([0,0])
  bb1=[]
  bb1.extend(df["base"])    
  bb1.extend([""])
  bb2=[]
  bb2.extend(df["base"])    
  bb2.extend(["",""])
  b=[]
  b.extend(df["base"]+bb1[1:]+bb2[2:])
  #b.extend(["",""])
  df["ref"]=b
  nu=[]
  nu.extend((df["P("+comb[0]+")"]*new1[2:]*new[1:]).round(3))
  nu.extend([0.000,0.00])
  li=[1.000,0.000,0.000]*(int(len(nu[2:])/3)+1)
  df[comb]=nu[:-2]
  df[comb]=df[comb]*li[0:len(df[comb])]
 els=[] 
 for line in mylist1.split("\n"):
    newd=line.split("=")[0].replace(" ","")[0]
    els.append(newd)
    ops=line.split("=")[1].split("+")
    df[newd]=df[ops[0].replace(" ","")]
    for el in ops[1:]:
      df[newd]= (df[newd]+ df[el.replace(" ","")]).round(3)
 df["newcol1"]=df[els].max(axis=1)
 df["alt"]=df[els].idxmax(axis="columns")
 tx=(df['pos']%3 ==r )
 emp=[]
 emp.extend(pi)
 emp.extend(els)
 emp.append("ref")
 emp.append("alt")
 emp.append("sorted")
 def add(arg):
    def floating(a):
       tt=args.treshold
       a0=int(a[0]*1000)/1000
       a1=a[1]
       return a0,a1  
    def get_proba(a):
        return a[0]
    arg=list(map(floating,arg))

    arg.sort(key=get_proba,reverse=True) 
    return arg
 def itere(row):
   return  add(list(zip(list(map(row.get,els)),els)))
 df['sorted'] = df.apply(itere , axis = 1)
 rr=len(df['sorted'])
 llo=list(range(2,(len(els)+3)*20000))
 df = df[cols]
 bahaee=[]
 for i in range(len(els)):
    try:   
     def getit(row):
       tt=args.treshold
       res=row.get("sorted")[i]
       #print(res[0],eval(tt))
       if res[0]>eval(tt):
           return res
       else:
          return ""

     bahae=" "*i
     df[bahae] = df.apply(getit , axis = 1)
     bahaee.append(bahae)

    except Exception as eror:
        print(eror)
 
 conteur=0
 
 verbos.extend(bahaee)
 try:
  if args.verbose:   
   df=df[verbos]
   print(111111111111111)
  else:
   newi=["gene","pos","ref","alt"]
   newi.extend(bahaee)
   df=df[newi]      
 except:
  print("looooooooooooooooool")   
  newi=["gene","pos","ref","alt"]
  newi.extend(bahaee)
  df=df[newi]
 for el,val in list(table.items()):
  df.loc[df['ref']==el.lower(),"ref"]=val
 
 #
 new=df[tx]
 new['pos'] = range(1,len(new['pos'])+1)
 
 new.to_csv(args.output_file_name, index=False)

