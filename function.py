
def returnSpeciseInfo():
    species = ['Acer_yangbiense', 'Amaranthus_hypochondriacus', 'Amborella_trichopoda', 'Aquilegia_coerulea', 'Arabidopsis_lyrata', 'Arabidopsis_thaliana', 'Arabis_nemorensis', 'Asparagus_officinalis', 'Beta_vulgaris', 'Boechera_stricta', 'Brachypodium_distachyon', 'Brassica_cretica', 'Brassica_napus', 'Brassica_oleracea', 'Brassica_rapa', 'Camelina_sativa', 'Camellia_sinensis', 'Cannabis_sativa', 'Capsella_grandiflora', 'Capsella_rubella', 'Capsicum_annuum', 'Capsicum_baccatum', 'Capsicum_chinense', 'Carica_papaya', 'Carpinus_fangiana', 'Chenopodium_quinoa', 'Cinnamomum_micranthum', 'Citrus_clementina', 'Citrus_sinensis', 'Citrus_unshiu', 'Coffea_arabica', 'Coffea_canephora', 'Coffea_eugenioides', 'Corchorus_capsularis', 'Daucus_carota', 'Durio_zibethinus', 'Eucalyptus_grandis', 'Eutrema_salsugineum', 'Fragaria_vesca', 'Ginkgo_biloba', 'Glycine_max', 'Gossypium_arboreum', 'Gossypium_barbadense', 'Gossypium_darwinii', 'Gossypium_hirsutum', 'Gossypium_raimondii', 'Helianthus_annuus', 'Herrania_umbratica', 'Hevea_brasiliensis', 'Ipomoea_nil', 'Ipomoea_triloba', 'Jatropha_curcas', 'Juglans_regia', 'Lactuca_sativa', 'Linum_usitatissimum', 'Macleaya_cordata', 'Malus_baccata', 'Malus_domestica', 'Manihot_esculenta', 'Microthlaspi_erraticum', 'Mimulus_guttatus', 'Morella_rubra', 'Morus_notabilis', 'Musa_acuminata', 'Nelumbo_nucifera', 'Nicotiana_attenuata', 'Nicotiana_sylvestris', 'Nicotiana_tabacum', 'Nicotiana_tomentosiformis', 'Nymphaea_colorata', 'Nymphaea_thermarum', 'Olea_europaea', 'Oryza_brachyantha', 'Oryza_sativa', 'Panicum_virgatum', 'Papaver_somniferum', 'Parasponia_andersonii', 'Phaseolus_vulgaris', 'Physcomitrella_patens', 'Pistacia_vera', 'Populus_alba', 'Populus_euphratica', 'Populus_trichocarpa', 'Prunus_avium', 'Prunus_dulcis', 'Prunus_mume', 'Prunus_persica', 'Punica_granatum', 'Pyrus_bretschneideri', 'Quercus_lobata', 'Quercus_suber', 'Raphanus_sativus', 'Rhamnella_rubrinervis', 'Ricinus_communis', 'Rosa_chinensis', 'Salix_brachista', 'Salix_purpurea', 'Setaria_italica', 'Solanum_chilense', 'Solanum_lycopersicum', 'Solanum_pennellii', 'Solanum_tuberosum', 'Sorghum_bicolor', 'Spinacia_oleracea', 'Tarenaya_hassleriana', 'Theobroma_cacao', 'Trema_orientale', 'Vitis_riparia', 'Vitis_vinifera', 'Zea_mays']
    #print(species)
    return species



def read_pfam_result(jsonfile):
    import json
    with open(jsonfile, "r") as f:
        dictionary = json.load(f)
    return dictionary

def find_gene_family(targetspecies,pfam):
    species = returnSpeciseInfo()
    result_protein = []
    species_target = {}
    pfam_list = pfam.split(",")
    if targetspecies in species:
        print("Reading json data ...")
        species_dictionary = read_pfam_result("E:\work\plantFamily\database\json\pfamscan_result.json")
        print("Finished reading! ")
        protein_info = species_dictionary[targetspecies]
        print("Now is searching %s" % targetspecies + "...")
        for protein in protein_info.keys():
            #pfamID = protein_info[protein]["clanID"]
            pfamID = protein_info[protein]["pfamID"]
            pfamID = [i[:7] for i in pfamID]
            no_exist = [False for a in pfam_list if a not in pfamID]
            if no_exist:
                pass
            else:
                result_protein.append(protein)
        result_protein = list(set(result_protein))
        print("Finished searching!")
        print(result_protein)
        return result_protein
    if targetspecies == "all":
        print("Reading json data ...")
        species_dictionary = read_pfam_result("E:\work\plantFamily\database\json\pfamscan_result.json")
        print("Finished reading! ")
        for key in species_dictionary.keys():
            print("Now is searching %s" % key + "...")
            if key not in species_target:
                species_target[key] = []
                for protein in species_dictionary[key].keys():
                    #pfamID = species_dictionary[key][protein]["clanID"]
                    pfamID = species_dictionary[key][protein]["pfamID"]
                    pfamID = [i[:7] for i in pfamID]
                    no_exist = [False for a in pfam_list if a not in pfamID]
                    if no_exist:
                        pass
                    else:
                        species_target[key].append(protein)
            else:
                for protein in species_dictionary[key].keys():
                    pfamID = species_dictionary[key][protein]["clanID"]
                    no_exist = [False for a in pfam_list if a not in pfamID]
                    if no_exist:
                        pass
                    else:
                        species_target[key].append(protein)
            species_target[key] = list(set(species_target[key]))
            print("Finished searching!")
        print(species_target)
        return species_target


def Multi_alignment(fastafile):
    from Bio.Align.Applications import MuscleCommandline
    muscle_exe = r"E:\work\plantFamily\soft\muscle3.8.31_i86win32.exe"
    inputfile = fastafile
    outfile = fastafile + ".aligned.aln"
    muscle_cline = MuscleCommandline(muscle_exe, input=inputfile, out=outfile, clwstrict=True)
    muscle_cline()

def find_gene_family_Seq(targetspecies,pfam,savepath):
    from Bio import Entrez
    Entrez.email = "A.N.Other@example.com"
    genefamily = find_gene_family(targetspecies,pfam)
    protein_info = {}
    result = open(savepath,"a")
    if isinstance(genefamily,list) == True:
        for i in genefamily:
            print("Searching %s protein sequence..." %i)
            try:
                protein_info[i] = ""
                record = Entrez.efetch(db="protein", id=i, rettype="fasta", retmode="text")
                fasta_txt = record.read()
                result.write(fasta_txt)
            except:
                errorfile = open(savepath + ".error.txt", "a")
                print(i + " maybe can find in Phytozome or NCBI, because this may not be annotated.")
                errorfile.write(i + "\n")
    else:
        for key in genefamily:
            if key not in protein_info.keys():
                print("Searching %s species..." % key)
                protein_info[key] = {}
                for protein in genefamily[key]:
                    try:
                        print("Searching %s protein sequence..." % protein)
                        protein_info[key][protein] = ""
                        record = Entrez.efetch(db="protein", id=protein, rettype="fasta", retmode="text")
                        fasta_txt = record.read()
                        result.write(fasta_txt)
                    except:
                        errorfile = open(savepath + ".error.txt", "a")
                        print(i + " maybe can find in Phytozome or NCBI, because this may not be annotated.")
                        errorfile.write(protein + "\n")
            else:
                print("Something is wrong!")
    print("Protein sequences had writen.")
    try:
        Multi_alignment(savepath)
    except:
        print("MUSCLE error, please submit the fasta file to https://www.ebi.ac.uk/Tools/msa/muscle/#")


def protein2gene(targetspecies,pfam,savepath):
    genefamily = find_gene_family(targetspecies,pfam)
    gene_info = {}
    from Bio import Entrez
    import re
    Entrez.email = "A.N.Other@example.com"
    if isinstance(genefamily,list) == True:
        for i in genefamily:
            print(i)
            print("Searching %s gene sequence..." %i)
            try:
                record = Entrez.efetch(db="protein", id=i, rettype="gb", retmode="text")
                fasta_txt = record.read()
                geneID = re.findall('coded_by=(.*)',fasta_txt)
                geneID = geneID[0].split(":")
                geneID = geneID[0]
                geneID = re.sub('"',"",geneID)
                record = Entrez.efetch(db="nucleotide", id=geneID, rettype="fasta", retmode="text")
                fasta_txt = record.read()
                result = open(savepath + ".geneSeq","a")
                result.write(fasta_txt)

            except:
                errorfile = open(savepath+".error.txt","a")
                errorfile.write(i + "\n")
    else:
        for key in genefamily:
            if key not in gene_info.keys():
                print("Searching %s species..." % key)
                gene_info[key] = {}
                for protein in genefamily[key]:
                    try:
                        print("Searching %s gene sequence..." % protein)
                        gene_info[key][protein] = ""
                        # record = Entrez.efetch(db="protein", id=protein, rettype="fasta", retmode="text")
                        # fasta_txt = record.read()
                        # geneID = re.findall('GeneID:(\d+)"', fasta_txt)
                        # record = Entrez.efetch(db="nucleotide", id=geneID[0], rettype="fasta", retmode="text")
                        # fasta_txt = record.read()


                        record = Entrez.efetch(db="protein", id=protein, rettype="gb", retmode="text")
                        fasta_txt = record.read()
                        geneID = re.findall('coded_by=(.*)', fasta_txt)
                        geneID = geneID[0].split(":")
                        geneID = geneID[0]
                        geneID = re.sub('"', "", geneID)
                        record = Entrez.efetch(db="nucleotide", id=geneID, rettype="fasta", retmode="text")
                        fasta_txt = record.read()
                        result = open(savepath, "a")
                        result.write(fasta_txt)

                    except:
                        errorfile = open(savepath + ".error.txt", "a")
                        print(key + "\t" + protein + " maybe can find in Phytozome or NCBI, because this may not be annotated.")
                        errorfile.write(key + "\t" + protein + "\n")

            else:
                print("Something is wrong!")
    try:
        Multi_alignment(savepath+".aligned.aln")
    except:
        print("MUSCLE error, please submit the fasta file to https://www.ebi.ac.uk/Tools/msa/muscle/#")