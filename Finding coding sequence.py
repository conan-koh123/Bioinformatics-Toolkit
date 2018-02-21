from Tkinter import *

root = Tk()
root.wm_title("Finding the coding sequence") #title of program
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

#Creating Protein label
protein_label = Label(top_frame, text = 'Protein of Interest:')
protein_label .grid(row=0, column=7, sticky=W, padx=4, pady = 4)

#Creating Protein textbox with scrollbar
protein_tb = Text(top_frame, height=7, width=35)
protein_tb.config(font=("Helvetica", 10), undo=True, wrap='word')
protein_tb.grid(row=1, column=7, columnspan = 5, rowspan = 5, padx=4, pady = 4)
scrollbar_protein=Scrollbar(top_frame, command=protein_tb.yview)
scrollbar_protein.grid(row=1, column= 13, rowspan=5, sticky=N+S+W)
protein_tb['yscrollcommand'] = scrollbar_protein.set

#Create protein name entry box
protein_name = Entry(top_frame)
protein_name.grid(row = 0, column = 8)
protein_name.insert(END, " ...(Name of protein) ...")

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
    (TAA, TAG, TGA) = len(DNA),len(DNA) ,len(DNA) 
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
    return protein_output[0:len(protein_output)]

#create a function to obtain coding sequence from protein of interest
def coding_sequence(DNA, protein):
    if protein[0] != "M":
        return "Problem: Protein of interest does not start with methionine"
    else:
        DNA_gene = DNA
        pos_codon_list = find_all_codon(DNA_gene)
        for pos_codon in pos_codon_list:
            DNA_gene = DNA_gene[pos_codon:]
            (TAA, TAG, TGA) = stop_codon(DNA_gene)
            min_stop_codon = min(TAA, TAG, TGA)
            coding_string = DNA_gene[0:min_stop_codon]
            protein_string = coding_to_protein(coding_string)
            if protein.replace("\n", ""). replace(" ", "") == protein_string:
                return DNA_gene[0:min_stop_codon + 3]
            else:
                DNA_gene = DNA_gene[3:]
        return "No coding sequence based on protein of interest"

#Create a function to print output in the output box
def output_box():
    output_tb.delete("1.0",END)
    try:
        DNA = gene_tb.get("1.0",END).replace("\n", "").replace(" ", "")
        protein = protein_tb.get("1.0",END).replace("\n", "").replace(" ", "")
        output_nucleotide = coding_sequence(DNA, protein)
        protein_output = str(protein_name.get()) + ": " + "\n" + protein + "\n"
        output_tb.insert(END, protein_output)
        seq_output =  "\n" + "Coding sequence" + "\n" + output_nucleotide + "\n"
        output_tb.insert(END, seq_output)
        if output_nucleotide != "No coding sequence based on protein of interest":
            before_coding = DNA[0:DNA.find(output_nucleotide)]
            before_coding_output = "\n" + str(gene_name.get()) + ":" + "\n" + before_coding
            output_tb.insert(END, before_coding_output)
            output_tb.insert(END, "{" + output_nucleotide + "}", 'Red')
            output_tb.tag_config('Red', foreground='red')
            after_coding = DNA[DNA.find(output_nucleotide) + len(output_nucleotide):]
            output_tb.insert(END, after_coding)
    except ValueError:
        output_error = "Error: Problem with one of the data fields"
        output_tb.insert(END, output_error)

#function allows selection of text by Ctrl + A for gene_tb
def select_all_gene_tb(event):
    gene_tb.tag_add(SEL, "1.0", END)
    gene_tb.mark_set(INSERT, "1.0")
    gene_tb.see(INSERT)
    return 'break'

#function allows selection of text by Ctrl + A for protein_tb
def select_all_protein_tb(event):
    protein_tb.tag_add(SEL, "1.0", END)
    protein_tb.mark_set(INSERT, "1.0")
    protein_tb.see(INSERT)
    return 'break'

#function allows selection of text by Ctrl + A for output_tb
def select_all_output_tb(event):
    output_tb.tag_add(SEL, "1.0", END)
    output_tb.mark_set(INSERT, "1.0")
    output_tb.see(INSERT)
    return 'break'
        
#Creating button to print output   
print_output = Button(top_frame, text="Print Output", command = output_box)
print_output .grid(row=1, column=16,sticky=E, padx = 5)    
    

gene_tb.bind("<Control-Key-a>", select_all_gene_tb) # bind event to select_all for Gene of interest textbox
protein_tb.bind("<Control-Key-a>", select_all_protein_tb) # bind event to select_all for protein of interest textbox
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
type instruction in output
reset button
try putting tabs
"""


