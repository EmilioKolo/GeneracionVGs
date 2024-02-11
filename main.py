from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq, MutableSeq

# Esto es por si se necesita cambiar el formato de aminoacidos
diccionario_aa = {"A":"Ala","C":"Cys","D":"Asp","E":"Glu","F":"Phe",
                  "G":"Gly","H":"His","I":"Ile","K":"Lys","L":"Leu",
                  "M":"Met","N":"Asn","P":"Pro","Q":"Gln","R":"Arg",
                  "S":"Ser","T":"Thr","V":"Val","W":"Trp","Y":"Tyr"}; 

# Primero se obtiene la secuencia de aminoacidos de referencia
Filename_prot_ref = 'P52952'; 
ext = 'fasta'; 
for seq_record in SeqIO.parse(Filename_prot_ref + '.' + ext, 'fasta'):
    NKX2_5_prot_ref = str(seq_record.seq);

# Despues se obtiene la secuencia de nucleotidos de referencia
Filename_exons = "NM_004387.4.exons.fa"; 

# NKX2_5_exons es la secuencia de los exones, incluyendo UTR
NKX2_5_exons = Seq("", IUPAC.unambiguous_dna);

for seq_record in SeqIO.parse(Filename_exons, "fasta"):
    NKX2_5_exons = NKX2_5_exons + str(seq_record.seq);

# Se eliminan a mano los UTR para tener la secuencia de referencia
NKX2_5 = NKX2_5_exons[123:-458];

# Se traduce para confirmar que es igual a la secuencia proteica de referencia
NKX2_5_prot = NKX2_5.translate(to_stop=True);

if str(NKX2_5_prot)==NKX2_5_prot_ref:
    print('La proteina traducida y la provista en ' + Filename_prot_ref + '.' + ext +
          ' son iguales.');

# Defino una funcion que devuelve todas las GVs que se pueden obtener con una sola variante
def devolverGVsSNP(ref):
    # Se le dan tres nucleotidos de referencia y prueba todas las SNP
    # Devuelve la lista de todos los cambios posibles
    GVs_SNP = [];
    aa_ref = ref.translate();
    # Un ciclo para cada nucleotido
    for nuc in range(len(ref)):
        # Un ciclo para cada posible variante del nucleotido
        for change in ["A","C","T","G"]:
            GV = MutableSeq(str(ref), IUPAC.unambiguous_dna);
            GV[nuc] = change;
            aa_GV = str(GV.toseq().translate(to_stop=True));
            # Prueba si la variante cambia el aminoacido
            if aa_GV!='' and aa_GV!=str(aa_ref) and not (aa_GV in GVs_SNP):
                GVs_SNP.append(aa_GV);
    return(GVs_SNP)


L = [];
for aa in range(len(NKX2_5_prot_ref)):
    # L es una matriz con una lista por posicion de AA
    # Cada posicion tiene una listas de posibles AAs a los que se llega con SNP
    ref = NKX2_5[aa*3:(aa+1)*3]; # Codon para el aminoacido
    L.append(devolverGVsSNP(ref));

# Version vieja de generarListaGVs
'''
lista_variantes = [];
HD_range = [-1,-1];
for aa in range(len(L)):
    pos = aa+1;

    if pos == HD[0]:
        HD_range[0] = len(lista_variantes);
    if pos == HD[1]+1:
        HD_range[1] = len(lista_variantes);

    ref_aa = NKX2_5_prot_ref[aa]; # Se puede usar diccionario_aa para codigo de 3 letras

    for unique_GV in L[aa]:
        new_aa = unique_GV; # Se puede usar diccionario_aa para codigo de 3 letras
        lista_variantes.append(str(ref_aa) + str(pos) + str(new_aa));

lista_HD = lista_variantes[HD_range[0]:HD_range[1]]
'''

# Funcion que transforma L en una lista de variantes en formato para FoldX
def generarListaGVs(cristal_range,GV_list,chain,prot_seq,threeletter=False,shift=0,mods=[]):
    diccionario_aa = {"A":"Ala","C":"Cys","D":"Asp","E":"Glu","F":"Phe",
                      "G":"Gly","H":"His","I":"Ile","K":"Lys","L":"Leu",
                      "M":"Met","N":"Asn","P":"Pro","Q":"Gln","R":"Arg",
                      "S":"Ser","T":"Thr","V":"Val","W":"Trp","Y":"Tyr"};
    L_GVs = GV_list;
    lista_variantes = [];
    list_range = [-1,-1];
    mod_p_seq = str(prot_seq);

    # Asigno las modificaciones a la secuencia proteica
    for i in mods:
        mod_pos = int(i[0]);
        mod_aa = str(i[1]);
        # Si la modificacion esta en L_GVs se hace un enroque simple
        if mod_aa in L_GVs[mod_pos-1]:
            mod_p_seq = mod_p_seq[:mod_pos-1] + mod_aa + mod_p_seq[mod_pos:];
            L_GVs[mod_pos-1][L_GVs[mod_pos-1].index(mod_aa)] = str(prot_seq)[mod_pos-1];
        # Si la modificacion no esta en L_GVs y no es el valor de referencia
        # Se agrega el valor de referencia a L_GVs
        elif mod_aa != mod_p_seq[mod_pos-1]:
            if not str(prot_seq)[mod_pos-1] in L_GVs[mod_pos-1]:
                L_GVs[mod_pos-1].append(str(prot_seq)[mod_pos-1]);
            mod_p_seq = mod_p_seq[:mod_pos-1] + mod_aa + mod_p_seq[mod_pos:];
        # Si la modificacion es igual al valor de referencia no se hace nada

    # Recorro la lista de variantes por posicion
    for aa in range(len(L_GVs)):
        pos = aa+1;
        # Defino el rango en lista_variantes
        if pos == cristal_range[0]:
            list_range[0] = len(lista_variantes);
        if pos == cristal_range[1]+1:
            list_range[1] = len(lista_variantes);
        if threeletter:
            for GV in L_GVs[aa]:
                lista_variantes.append(diccionario_aa[str(mod_p_seq[aa])] + str(chain) +
                                       str(pos+shift) + diccionario_aa[str(GV)] + ';');
        else:
            for GV in L_GVs[aa]:
                lista_variantes.append(str(mod_p_seq[aa]) + str(chain) +
                                       str(pos+shift) + str(GV) + ';');
    R = lista_variantes[list_range[0]:list_range[1]];
    return(R)


# Una vez hecha la funcion, armo listas de valores para correr
# Para cada cristal, tomo una cadena, rango de residuos, desfasaje y modificaciones
lista_rangos_cristales = [('3rkq',
                           [(137,194,'A',0,[(137,'G'),(193,'S')]),
                            (137,193,'B',0,[(137,'G'),(193,'S')])]),
                          ('4s0h',
                           [(142,194,'B',0,[(193,'S')]),
                            (142,194,'F',-37,[(193,'S')])]),
                          ('Modeller',
                           [(137,194,'A',-136,[])])
                          ];

# Guardo los resultados de la funcion en un archivo por cristal
for cristal in lista_rangos_cristales:
    F = open('../variables/individual_list' + cristal[0] + '.txt','w');
    F.close();

    for cad in range(len(cristal[1])):
        cadena = cristal[1][cad];
        lista_HD = generarListaGVs((cadena[0],cadena[1]),L,cadena[2],NKX2_5_prot_ref,
                                   shift=cadena[3],mods=cadena[4])
        with open('../variables/individual_list' + cristal[0] + '.txt','a') as archivo_texto:
            for i in range(len(lista_HD)):
                if i==len(lista_HD)-1 and cad==len(cristal[1])-1:
                    archivo_texto.write(lista_HD[i]);
                else:
                    archivo_texto.write(lista_HD[i]);
                    archivo_texto.write("\n");

print('Programa finalizado. Revisar la carpeta "variables" para ver los archivos.')
