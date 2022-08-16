from Bio import SeqIO
from Bio import Entrez
from Bio import ExPASy
from Bio.Seq import Seq  ## IUPAC alphabet
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord

path = str(
    "C:\\Users\putnt\Documents\Arin's Folder\BACTERIOPHAGE GENOMIC DATA\Bacteriophage contigs\Redo\SESP1-500.fas")

for seq_record in SeqIO.parse(path, "fasta", alphabet=IUPAC. ambiguous_dna):  # foreach loop
    print("id: " + seq_record.id)
    print("repr: " + repr(seq_record.seq))
    print("length: ", (len(seq_record)))

length, id, sequence for each contig
converts single alpha to IUPAC

### Dictionary keys
# sepsDic= SeqIO.to_dict(SeqIO.parse(path, "fasta", alphabet=IUPAC.ambiguous_dna))
# path.close()
# print(sepsDic.keys())
# ##Trying to specify dictionary keys

all_species = []
for seq_record in SeqIO.parse(path, "fasta") :
    all_species.append(seq_record.description.split()[1])
path.close()
print (all_species)

###Converting fasta file to genbank 4.4.1 (reverse)
inputFasta = open(path, "rU")
outputGenbank = open("seps1.gbk", "w")

sequences = list(SeqIO.parse(inputFasta, "fasta"))

for seq in sequences:
    seq.seq.alphabet = generic_dna

count = SeqIO.write(sequences, outputGenbank, "genbank")

# outputGenbank.close()
inputFasta.close()
print("Converted %i records" % count)
#
# outputGenbank.read()


###Reading  gebank file
seps1gb = open("seps1.gbk")
for seq_record in SeqIO.parse(seps1gb, "genbank"):
    print("id: " + seq_record.id)
    print("repr: " + repr(seq_record.seq))
    print("length: ", (len(seq_record)))

###Converting a file of seq to comp.
### extracting data
#
# recordIterator=SeqIO.parse(seps1gb, "genbank")
# firstRecord= recordIterator.next()
# print(firstRecord)

###Reverse complement conversion
for record in SeqIO.parse(seps1gb, "genbank"):
    print(record.id)
    print(record.seq.reverse_complement().tostring())
seps1gb.close()

##save reverse comp file
def make_rc_record(record):
    """Return a new SeqRecord with the reverse complement sequence."""
    rc_rec = SeqRecord(seq=record.seq.reverse_complement(), \
                       id="rc_" + record.id, \
                       name="rc_" + record.name, \
                       description="reverse complement")
    return rc_rec
