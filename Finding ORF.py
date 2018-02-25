from Tkinter import *
import tkFileDialog
import tkMessageBox

root = Tk()
root.wm_title("Finding Open Reading Frames (ORF)") #title of program
root.geometry("700x700") #size of root (canvas)

#Create top and bottome frame
top_frame = Frame(root, bg='lightblue', width=600, height=400, padx=3)
btm_frame = Frame(root, width=600, height=300, padx=3)

root.grid_rowconfigure(1, weight=1)
root.grid_columnconfigure(0, weight=1)

top_frame.grid(row=0, sticky="ew")
btm_frame.grid(row=1, sticky="ew")

#Creating Gene label
gene_label = Label(top_frame, text = 'mRNA of Interest:')
gene_label .grid(row=0, column=0, sticky=W, padx=4, pady = 4)

#Creating Gene textbox with scrollbar
gene_tb = Text(top_frame, height=7, width=35)
gene_tb.config(font=("Helvetica", 10), undo=True, wrap='word')
gene_tb.grid(row=1, column=0, columnspan = 5, rowspan = 5, padx=4, pady = 4)
scrollbar_gene=Scrollbar(top_frame, command=gene_tb.yview)
scrollbar_gene.grid(row=1, column=6, rowspan=5, sticky=N+S+W)
gene_tb['yscrollcommand'] = scrollbar_gene.set

#Create Gene name entry box
gene_name = Entry(top_frame)
gene_name.grid(row = 0, column = 1)
gene_name.insert(END, " ...(Name of mRNA) ...")

#Creating Output label
output_label = Label(btm_frame, text = 'Output:')
output_label .grid(row=0, column=0, sticky=W, padx=4, pady = 4)

#Creating Output textbox with scrollbar
output_tb = Text(btm_frame, height=30, width=93)
output_tb.config(font=("Helvetica", 10), undo=True, wrap='word')
output_tb.grid(row=1, column=0, columnspan = 20, rowspan = 10,
               padx=4, pady = 4, sticky = W)
scrollbar_output=Scrollbar(btm_frame, command=output_tb.yview)
scrollbar_output.grid(row=1, column= 21, rowspan=10, sticky=N+S+W)
output_tb['yscrollcommand'] = scrollbar_output.set

#Insert output textbox with instruction
instruction = "Instructions to using 'Finding Open Reading Frame (ORF)' program: \n \
1) Enter mRNA of interest and name of mRNA. Ensure that \
non 'ATCG' characters are not entered. \n \
2) Press 'Print Output' for the result. Section 1 contains all ORFs from \
mRNA, while Section 2 contains all ORFs from its reverse complement \n \
3) Press 'Save Output' to save the result as a text file."

output_tb.insert(END, instruction)

#protein_dictionary
protein_dictionary = {"TTT": "F",
                         "TTC": "F", 
                         "TTA": "L", 
                         "TTG": "L", 
                         "CTT": "L", 
                         "CTC": "L",
                         "CTA": "L",
                         "CTG": "L",
                         "ATT": "I", 
                         "ATC": "I",
                         "ATA": "I", 
                         "ATG": "M", 
                         "GTT": "V",
                         "GTC": "V", 
                         "GTA": "V", 
                         "GTG": "V",
                         "TCT": "S", 
                         "TCC": "S",
                         "TCA": "S",
                         "TCG": "S", 
                         "CCT": "P", 
                         "CCC": "P", 
                         "CCA": "P", 
                         "CCG": "P", 
                         "ACT": "T", 
                         "ACC": "T", 
                         "ACA": "T", 
                         "ACG": "T", 
                         "GCT": "A",
                         "GCC": "A", 
                         "GCA": "A", 
                         "GCG": "A", 
                         "TAT": "Y", 
                         "TAC": "Y", 
                         "TAA": "-", 
                         "TAG": "-", 
                         "CAT": "H",
                         "CAC": "H", 
                         "CAA": "Q", 
                         "CAG": "Q", 
                         "AAT": "N", 
                         "AAC": "N", 
                         "AAA": "K", 
                         "AAG": "K", 
                         "GAT": "D",
                         "GAC": "D",
                         "GAA": "E", 
                         "GAG": "E", 
                         "TGT": "C", 
                         "TGC": "C", 
                         "TGA": "-",
                         "TGG": "W", 
                         "CGT": "R", 
                         "CGC": "R", 
                         "CGA": "R", 
                         "CGG": "R", 
                         "AGT": "S", 
                         "AGC": "S", 
                         "AGA": "R", 
                         "AGG": "R", 
                         "GGT": "G", 
                         "GGC": "G",
                         "GGA": "G", 
                         "GGG": "G"}

#create a function to find all codon position 
def find_all_codon(gene):
    pos_list = []
    pos_gene = gene.find("ATG")
    while pos_gene != -1:
        pos_list.append(pos_gene)
        gene = gene[pos_gene+3:]
        pos_gene = gene.find("ATG")
    return pos_list

#create a function to find stop codon in multiples of 3
def stop_codon(DNA):
    (TAA, TAG, TGA) = len(DNA) + 2,len(DNA) + 2 ,len(DNA) + 2 
    for index in range(0, len(DNA), 3):
        if DNA[index: index + 3] == "TAA":
            TAA = index
            break
        elif DNA[index: index + 3] == "TAG":
            TAG = index
            break
        elif DNA[index: index + 3] == "TGA":
            TGA = index
            break
    return (TAA, TAG, TGA)

#create a function for translating coding sequence into protein
def coding_to_protein(gene):
    starting_position = gene.find("ATG")
    gene = gene.replace("/n", "")
    translating_seq = gene[starting_position:]
    translating_sequence_1 = translating_seq[0: len(translating_seq) // 3 * 3]
    protein_output = ""
    for nucleotide in range(0, len(translating_sequence_1), 3):
        amino_acid_op =  protein_dictionary[translating_sequence_1[nucleotide:nucleotide + 3]]
        if amino_acid_op != "Stop":
            protein_output = protein_output + amino_acid_op 
        else: 
            break
    return protein_output[0:len(protein_output) -1]

#create a function to reverse complement a mRNA
def mrna_reverse_complement(mrna):
    mrna_complement = ""
    for nucleotide in mrna:
        if nucleotide == "A":
            mrna_complement = mrna_complement + "T"
        elif nucleotide == "T":
            mrna_complement = mrna_complement + "A"
        elif nucleotide == "C":
            mrna_complement = mrna_complement + "G"
        elif nucleotide == "G":
            mrna_complement = mrna_complement + "C"
    return mrna_complement[::-1]

#Create a list of ORF
def ORF(DNA):
    orf_list = []
    pos_atg = DNA.find("ATG")
    while pos_atg != -1:
        DNA = DNA[pos_atg:]
        (TAA, TAG, TGA) = stop_codon(DNA)
        min_stop_codon = min(TAA, TAG, TGA)
        if min_stop_codon != len(DNA) + 2:
            sequence = DNA[0:min_stop_codon + 3]
            orf_list.append(sequence)
        #print("Sequence to list:" + sequence)
        DNA = DNA[3:]
        #print("sequence after:" + DNA)
        pos_atg = DNA.find("ATG")
        #print("pos_atg:" + str(pos_atg))
    return orf_list

#Create a function to insert ORF and protein sequence for mRNA and reverse complement
def printing_function():
    output_tb.delete("1.0",END)
    orf_original_sequence = gene_tb.get("1.0",END).replace("\n", "").replace(" ", "")
    list_orf_original = ORF(orf_original_sequence)
    section_1 = "\n" + "///////////////////////////////////////////////////////////////////////////////////" + "\n" + \
                "SECTION 1: ORF from " + gene_name.get() + "\n" + \
                "///////////////////////////////////////////////////////////////////////////////////" +\
                "\n"
    output_tb.insert(END, section_1)
    if list_orf_original != []:
        for orf in list_orf_original:
            orf_index = orf_original_sequence.find(orf)
            orf_insert = "\n" +  "ORF " + str(list_orf_original.index(orf) + 1) + ":" +  \
                         "\n" + "From position " + str(orf_index + 1) + " - " + str(orf_index + len(orf)) + \
                         "\n" + \
                         orf
            output_tb.insert(END, orf_insert)
            protein_original = coding_to_protein(orf)
            protein_insert = "\n" + "Protein " + str(list_orf_original.index(orf) + 1)  + ":" + \
                             "\n" + protein_original + "\n"
            output_tb.insert(END, protein_insert)
            output_tb.insert(END, " ==========================")
    else:
        output_tb.insert(END, "No ORFs are found for mRNA")
    section_2 =  "\n" + "///////////////////////////////////////////////////////////////////////////////////" + "\n" + \
                "SECTION 2: ORF from reverse complement " + gene_name.get() + "\n" + \
                "///////////////////////////////////////////////////////////////////////////////////" +\
                "\n"
    output_tb.insert(END, section_2)
    orf_reverse_complement = mrna_reverse_complement(orf_original_sequence)
    list_orf_reverse = ORF(orf_reverse_complement)
    count = 1
    if list_orf_reverse != []:
        for orf in list_orf_reverse:
            orf_index = orf_reverse_complement.find(orf)
            orf_insert = "\n" +  "ORF " +  str(count) + " :" + "\n" + \
                         "From position " + str(orf_index + 1) + " - " + str(orf_index + len(orf)) + \
                         "\n" + \
                         orf
            output_tb.insert(END, orf_insert)
            protein_reverse = coding_to_protein(orf)
            protein_insert = "\n" + "Protein " + str(count)  + ":" + \
                             "\n" + protein_reverse  + "\n"
            output_tb.insert(END, protein_insert)
            output_tb.insert(END, " ==========================")
            count += 1
    else:
        output_tb.insert(END, "No ORFs are found for reverse complement")


#function allows selection of text by Ctrl + A for gene_tb
def select_all_gene_tb(event):
    gene_tb.tag_add(SEL, "1.0", END)
    gene_tb.mark_set(INSERT, "1.0")
    gene_tb.see(INSERT)
    return 'break'


#function allows selection of text by Ctrl + A for output_tb
def select_all_output_tb(event):
    output_tb.tag_add(SEL, "1.0", END)
    output_tb.mark_set(INSERT, "1.0")
    output_tb.see(INSERT)
    return 'break'

#function to reset everything
def reset():
    output_tb.delete("1.0",END)
    gene_tb.delete("1.0",END)
    gene_name.delete(0,END)
    protein_name.delete(0,END)
    gene_name.insert(END, " ...(Name of mRNA) ...")
    protein_name.insert(END, " ...(Name of protein) ...")
        
#Creating button to print output   
print_output = Button(top_frame, text="Print Output", command = printing_function)
print_output .grid(row=1, column=16,sticky=E, padx = 5)


#Creating button to reset   
reset = Button(top_frame, text="Reset", command = reset)
reset.grid(row=3, column=16,sticky=W, padx = 5)    
    

#Create a function to save From:http://knowpapa.com/text-editor/
def save_command():
    file = tkFileDialog.asksaveasfile(mode='w', defaultextension=".txt")
    if file != None:
    # slice off the last character from get, as an extra return is added
        data = output_tb.get('1.0', END+'-1c')
        file.write(data)
        file.close()


#Creating button to save  
save = Button(top_frame, text="Save Output", command = save_command)
save.grid(row=1, column=19,sticky=W, padx = 5)    

gene_tb.bind("<Control-Key-a>", select_all_gene_tb) # bind event to select_all for Gene of interest textbox
output_tb.bind("<Control-Key-a>", select_all_output_tb) # bind event to select_all for output textbox

root.mainloop()


"""
Gene used for testing
1)TPH 1 (NM_004179.2)
2) Slc2a2 (NM_031197.2)
3) ESR1 (NM_00125.3)
4) GHR mRNA (NM_00163.4)
5) MCHR1 (NM_005297.3)
6) GRIA1 (NM_000827.3)
to do:
try putting tabs
"""
