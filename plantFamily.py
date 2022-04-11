import optparse
import  function
from optparse import OptionParser
from threading import Thread

# 启动子由于Selaginella_moellendorffii和Physcomitrella_patens基因组组装不完全，只有基因组草图，所以无法提取启动子序列！
def main():
    usage = "Usage: python plantFamily.py [options] arg1 arg2 ..."
    parser = optparse.OptionParser(usage, version="plantFamily 1.0")
    # 1.返回物种信息
    parser.add_option('--species-info', dest='speciesInfo', action="store_true",
                      help='Return species info in the database.')
    # 2.返回基因家族
    parser.add_option("--genefamily", dest="genefamily", action='store_true',
                      help="If you specify a species, it will return the gene family of that species. Otherwise, it will return the gene family of 109 species.")
    parser.add_option('-s','--species', dest="species", type='string', help="Your target species. If 'all' is chosed, it will return the gene family of 109 species.")
    parser.add_option('-p','--pfam', dest='pfam', type='string',help='Your query protein domain CLAN    .')
    # 3.返回查询蛋白的同源蛋白序列
    #parser.add_option("--genefamilySeq", dest="genefamilySeq", action='store_true',help="Returns the sequence of the homologous protein.")
    parser.add_option("--homologousSeq", dest="homologousSeq", action='store_true',help="Returns the sequence of the homologous protein and give a alignment result by muscle. -S is required.")
    parser.add_option('-S','--savepath', dest='savepath', type='string',help='Save path.')
    # 4.查找蛋白对应的基因
    parser.add_option('-g', "--protein2gene", dest="protein2gene", action='store_true',
                      help="Returns gene name and gene sequence.")


    (options, args) = parser.parse_args()

    # 1.返回物种信息的主程序
    if options.speciesInfo == True:
        species = function.returnSpeciseInfo()
        print(species)
    # 2.返回基因家族信息, -p PFXXXXXX
    if options.genefamily == True:
        parser = optparse.OptionParser("usage: python plantFamily.py --genefamily " + "-s <species> -p <pfamID,pfamID,pfamID...> [option]")
        if (options.species==None) | (options.pfam==None):
            print(parser.usage)
        # if options.savepath != None:
        #     targetspecies = options.species
        #     pfam = options.pfam
        #     savepath = options.savepath
        #     result_species_homologProtein = function.find_gene_family(targetspecies, pfam)
        #     result = open(savepath,"a")
        #     result.write("species\tproteinID\n")
        #     if options.species == "all":
        #         for key in result_species_homologProtein:
        #             for id in result_species_homologProtein[key]:
        #                 result.write(key + "\t" + id + "\n")
        if options.homologousSeq != None:
            targetspecies = options.species
            pfam = options.pfam
            savepath = options.savepath
            function.find_gene_family_Seq(targetspecies, pfam, savepath)
        if options.protein2gene != None:
            targetspecies = options.species
            pfam = options.pfam
            savepath = options.savepath
            function.protein2gene(targetspecies, pfam, savepath)

        if (options.homologousSeq == None) and (options.protein2gene == None) :
            targetspecies = options.species
            pfam = options.pfam
            t = Thread(target=function.find_gene_family, args=(targetspecies, pfam))
            t.start()

    # 3.返回同源蛋白序列 (在2.0版本中已经删去)
    # if options.genefamilySeq == True:
    #     parser = optparse.OptionParser("usage: python plantFamily.py --genefamilySeq " + "-s <species> -p <pfamID,pfamID,pfamID...> -S <save path>")
    #     if (options.species==None) | (options.pfam==None) | (options.savepath==None):
    #         print(parser.usage)
    #     else:
    #         targetspecies = options.species
    #         pfam = options.pfam
    #         savepath = options.savepath
    #         t = Thread(target=function.find_gene_family_Seq, args=(targetspecies, pfam, savepath))
    #         t.start()

    # 4.查找蛋白对应的基因
    # if options.protein2gene == True:
    #     parser = optparse.OptionParser("usage: python plantFamily.py --protein2gene " + "-s <species> -p <pfamID,pfamID,pfamID...> -S <save path>")
    #     if (options.species==None) | (options.pfam==None) :
    #         print(parser.usage)
    #     else:
    #         targetspecies = options.species
    #         pfam = options.pfam
    #         savepath = options.savepath
    #         t = Thread(target=function.protein2gene, args=(targetspecies, pfam, savepath))
    #         t.start()


if __name__ == "__main__":
    main()




