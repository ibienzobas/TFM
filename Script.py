#importa las librerias necesarias
import os
import getopt
import sys
from datetime import datetime



def get_transcripts_id(compound, dir): # Dado un archivo procedente del resultado de un tBLASTn, devuelve los identificadores de los transcritos que se encuentran en él.


    os.chdir(dir)

    genes_id = [] # Lista donde se guardan los identificadores de los transcritos

    for line in open(compound, "r", errors="ignore"): # Itera por cada linea (=transcrito) del fichero resultante del tBLASTn
        line = line.split()
        transcript = line[1].split(".") # Guardamos el identificador
        gene_id = ".".join(transcript[:2])
        genes_id.append(gene_id) # Lo introducimos en la lista anteriormente creada

    return(genes_id)



def get_genes_info(genes_id, compound_name, output_dir,genome): # Dado un conjunto de identificadores de transcritos, crea un archivo con la información de los genes que los codifican.

    os.chdir(output_dir)
    output_file = compound_name + ".genes" # Genera el archivo donde se guardan los posibles genes codificantes para los transcritos.
    with open(output_file, "a+") as f:
        for gene_id in genes_id: # Itera por cada identificador de transcrito
            for line in open(genome,"r", errors="ignore"): # Itera por cada linea del genoma anotado
                if line[0] != "#":
                    line=line.split()
                    if line[2] == "gene":
                        if line[8][3:] == gene_id:

                            f.write("   ".join(line) + "\n") # Si la linea corresponde con un gen, se guarda en el fichero ".genes"

    f.close()

    output_ruta=os.path.join(output_dir,output_file)

    return(output_ruta) # Devuelve la ruta con la localizacion del archivo ".genes", donde se ha guardado los posibles genes codificantes para los transcritos del input.



def get_mutations(genes_file,variants_dir,output_dir): # Dado un conjunto de genes, filtra las mutaciones que les afectan, en cada una de las variedades de melón (8).

    os.chdir(variants_dir)

    for variant in os.listdir(variants_dir):# Itera por cada archivo ".vcf" (8; uno de cada variedad de melón)

        if variant.split(".")[-1] == "vcf": # Comprueba que el archivo que va a utilizar es un vcf.
            variant_filename_path = os.path.join(output_dir, variant.split(".")[0]) # Crea el archivo que contiene las mutaciones para la variedad de melon iterada.
            with open(variant_filename_path + ".mutations","a+") as f: # Crea el fichero donde se guardan las mutaciones "on target" (en los genes de interés), para la variedad de melon iterada.
                os.chdir(output_dir)
                for line in open(genes_file, "r", errors="ignore"): # Itera por cada uno de los posibles genes codificantes de la enzima introducida.

                        line = line.split()
                        gene_id = line[8][3:] # Guarda el id del gen.
                        gene_pos=[line[0],line[3],line[4]] # Guarda la posicion del gen (cromosoma, posicion de inicio, posicion de fin).


                        f.write("\n" + "\n")

                        f.write(gene_id + ": mutations found in variant" + variant.split(".")[0] + "\n" + "\n") # Escribe en el archivo las mutaciones que vas a aencontrar para el gen iterado.

                        if len(gene_pos) != 0: # Si el archivo genes no esta vacio, continua.
                            variant_dir = os.path.join(variants_dir,variant)

                            for line_v in open(variant_dir, "r", errors="ignore"): # Abre el archivo vcf correspondiente.
                                if not line[0] == "#": # Omite las lineas que no contienen mutaciones (son comentarios).
                                    line_v=line_v.split()

                                    if line_v[0] == gene_pos[0]: # Si el cromosoma coincide
                                        if line_v[1] >= gene_pos[1] and line_v[1] <= gene_pos[2]: # Y la posicion de la mutacion esta entre la posicion de inicio y fin del gen iterado.
                                            f.write("\t".join(line_v) + "\n") # Guarda la mutacion en el archivo de salida.
                os.chdir(variants_dir)

            f.close()


def unique_mutations(compounds_dir, compound, variant): # Dado una enzima y una variedad de melon, filtra las mutaciones que afectan a los posibles genes de sintesis de dicha enzima. Deben ser unicas de esa variedad, es decir, no aparecen en las demas variedades de melon.


    wd=os.path.join(compounds_dir,compound)
    os.chdir(wd)
    all_mutations = [] # Lista donde se guarda las mutaciones que aparecen en la variedad de melon seleccioanda, en los posibles genes de sintesis de la enzima seleccionada.

    mutations_progresion={} # Destinado a la elaboracion de tablas que, para cada gen (gen_id), muestra el numero de mutaciones presentes en la variedad de melon para dicho gen, asi como las unicas de la variedad de melon seleccionada.

    for line in open(variant,"r"):
        line=line.split()
        if not len(line) == 0: # Omite las lineas vacias.

            if line[0][:4] == "MELO": # Para cada gen, guarda determinada informacion respecto a las mutaciones que le afectan.
                gen_id=line[0][:-1]
                mutations_progresion[gen_id]=[0,0,0,0]

            if line[0][:3] == "chr": # Si la linea es una mutacion, continua.

                mutations_progresion[gen_id][0] += 1 # Suma 1 al contador de mutaciones del gen correspondiente.
                all_mutations.append([line[0],line[1],line[3],line[4]]) # Guarda en la lista la informacion de esa mutacion (cromosoma, posicion, variante de referencia, alteracion).

    print("Total mutations found in ", compound ,"synthesis genes: ",len(all_mutations))  # Pon en pantalla el numero de mutaciones que afectan a la enzima y variedad de melon seleccionadas.

    for mutation_file in os.listdir(wd): # Itera por los archivos ".mutations" a excepcion del correspondiente a la variedad de melon seleccionada.

        if mutation_file.split(".")[-1] == "mutations": # comprueba que es un archivo ".mutations".

            if mutation_file.split(".")[0] != variant.split(".")[0]: #evita que itere sobre el archivo ".mutations" de la variante seleccionada, para la que buscamos mutaciones unicas.

                for line in open(mutation_file,"r"): # Abre el archivo ".mutations".
                    line=line.split()
                    if not len(line) == 0: # Omite las lineas vacias
                        if line[0][:3] == "chr": # Si la linea corresponde a una mutacion, continua.

                            posible_unique_mutation = [line[0],line[1],line[3],line[4]] # Guarda la informacion de la mutacion potencialmente unica.


                            if posible_unique_mutation in all_mutations: # Si encuentra la mutacion en la lista (que guarda las mutacuones de la variedad seleccionada), eliminala de la lista, dado que no es unica de la variedad de melon seleccionada.
                                all_mutations.remove(posible_unique_mutation)

    unique_mutations_in_variant = len(all_mutations) # Mutaciones restantes, que no han sido eliminadas, por lo que son unicas de la variedad de melon seleccionada.

    with open(variant.split(".")[0] + ".UniqueMutations", "w+") as f: # Crea el archivo donde vas a guardar las mutaciones unicas.

        for line in open(variant, "r"): # Itera por el archivo que contiene, para la enzima y variedad de melon seleccionadas, las mutaciones existentes.
            line = line.split()
            if not len(line) == 0: # Omite las lineas vacias.

                if line[0][:4] == "MELO": # Si la linea es la cabecera de un gen, haz lo siguiente.
                    gen_id=line[0][:-1] # Actualiza el valor del id del gen.
                    f.write("\n" + "\n" + "\t".join(line) + "\n" + "\n") # Escribe la linea en el archivo "unique_mutations".
                if line[0][:3] == "chr": # Si la linea es una mutacion

                    posible_unique_mutation= [line[0],line[1],line[3],line[4]] # Guarda la informacion referida a dicha mutacion.
                    if posible_unique_mutation in all_mutations: # Si la mutacion es unica, escribela en el archivo "unique_mutations".
                        mutations_progresion[gen_id][1] += 1 # Suma 1 al contador de mutaciones unicas para el gen correspondiente.
                        f.write("\t".join(line) + "\n")
    f.close()
    print(mutations_progresion)
    with open(variant.split(".")[0] + ".UniqueMutations", "a") as f: # Adjunta al final del archivo "unique_mutations" la tabla con los contadores de mutaciones totales y unicas para la enzima y variedad de melon seleccionadas.
        for gene, value in mutations_progresion.items():
            #f.write( gene + ": " + str(value[0]) + " mutaciones encontradas en la variante " + variant.split(".")[0] + ", " + str(value[1]) + " únicas." + "\n")
            f.write(gene + "\t" + str(value[0]) + "\t" + str(value[1]) + "\n")
    f.close()

    print("Unique mutations found for ", variant.split(".")[0],": ",len(all_mutations)) # Pon en pantalla el numero de mutaciones que afectan a la enzima y unicamente variedad de melon seleccionadas.





def _main_(argv): # Todo lo expuesto a paprtir de aqui esta unicamente relacionado con la forma en la que se lanza el script.

    if not argv:
        print("Please use the following input: ", "\n")
        print('script.py -b < BLASTs directory> -v <variants directory> -g <genome file> -o <output directory>')
        print('script.py -m <unique mutations finder mode> -d <compounds directory> -c <compound> -v <variant>')

        sys.exit(1)
    try:
      opts, args = getopt.getopt(argv,'hmb:v:o:g:d:')

    except getopt.GetoptError:
        print("Input format not valid")
        print('script.py -b < BLASTs directory> -v <variants directory> -g <genome file> -o <output directory>')
        print('script.py -m <unique mutations finder mode> -d <compounds directory> -c <compound> -v <variant>')

        sys.exit(2)

    BLASTresults_dir, variants_dir, genome, output_dir  = '', '', '', ''

    for opt, args in opts:
        if opt == '-m':

            compounds_dir = ''

            for opt, arg in opts:
                if opt == '-d':
                    compounds_dir = arg

            if compounds_dir == '':
                print("<<ERROR>> Compounds directory is missing")
                sys.exit(1)

            os.chdir(compounds_dir)
            compounds = {}
            c = 1

            for compound in os.listdir(compounds_dir):
                compounds[c] = compound
                c += 1

            for number, compound in compounds.items():
                print(number, compound, "\n")


            choice = int(input("Please choose a compound: "))

            try:
                compound = compounds[choice]
            except:
                print("Please choose a compound")
                sys.exit()

            #variants = {1: "Am_anot.mutations", 2: "Ch_anot.mutations", 3: "Gal_anot.mutations", 4: "PSA_anot.mutations", 5: "PSU_anot.mutations", 6: "Roch_anot.mutations", 7: "TB_anot.mutations", 8: "TN_anot.mutations" }
            #variants = {1: "Am.mutations", 2: "Ch.mutations", 3: "Gal.mutations", 4: "PSA.mutations", 5: "PSU.mutations", 6: "Roch.mutations", 7: "TB.mutations", 8: "TN.mutations" }
            variants = {1: "variants-Am.mutations", 2: "Ch.mutations", 3: "Gal.mutations", 4: "variants-PSA.mutations", 5: "PSU.mutations", 6: "Roch.mutations", 7: "variants-TB.mutations", 8: "TN.mutations" }


            for number, variant in variants.items():
                print(number, variant, "\n")
            choice = int(input("Please choose a variant: "))

            try:
                variant = variants[choice]
            except:
                print("Please choose a variant")
                sys.exit()


            unique_mutations(compounds_dir,compound,variant)
            sys.exit()

    for opt, arg in opts:
        if opt == '-h':
            print('script.py  -b < BLASTs directory> -v <variants directory> -g <genome file> -o <output directory>')
            print('script.py -m <unique mutations finder mode> -d <compounds directory> -c <compound> -v <variant>')
            sys.exit()
        if opt == '-b':
            BLASTresults_dir = arg
        elif opt == '-v':
            variants_dir = arg
        elif opt == '-g':
            genome = arg
        elif opt == '-o':
            output_dir = arg


    if BLASTresults_dir == '' or variants_dir == '' or genome == '' or output_dir == '':
        print("<<ERROR>> some input files are missing")
        sys.exit(1)

    print("Process starts... " + str(datetime.now()))
    os.chdir(BLASTresults_dir)

    for compound in os.listdir(BLASTresults_dir): # Para cada una de las enzimas presentes en el directorio de trabajo, haz lo siguiente.
        if compound[0] != ".": # Omite los archivos que comienzen por un punto.
            output_dir_1 = output_dir # Variable auxiliar
            output_dir_1 = os.path.join(output_dir_1, compound.split(".")[0]) # Crea una ruta con el directorio de trabajo + una carpeta con el nombre de la enzima iterada.
            os.mkdir(output_dir_1) # Crea dicho directorio.
            genes_id = get_transcripts_id(compound,BLASTresults_dir) # Guarda los genes que posiblemente codifican para dicha enzima.
            genes_file = get_genes_info(genes_id, compound, output_dir_1,genome) # Crea el archivo con la informacion de dichos genes.
            get_mutations(genes_file, variants_dir, output_dir_1) # Filtra las mutaciones que afectan a dichos genes, en cada variedad de melon (8 variedades de melon; por lo tanto crea 8 archivos distintos).



    print("Process finished... " + str(datetime.now()))

    sys.exit(0)

_main_(sys.argv[1:])