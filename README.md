# Trabajo-final-de-bases-ecol-gicas-y-gen-micas-de-interacci-n-organismo-ambiente

#Pasos #Acceso al servidor del lab de genética
ssh celiabc@132.248.15.41 -p 1967 -o ServerAliveInterval=60

#######Copiar las secuencias crudas desde la compu al servidor#####
1. Subir archivos .fastq.gz al servidor

En mi computadora, agrupar todos los datos crudos en una carpeta denominada BEGIOAM
Crear previamente en el servidor la carpeta BEGIOAM Para crear carpetas usar mkdir Para eliminar carpetas rmdir
Crear dentro de esta carpeta otra carpeta con el nombre Secuencias_crudas (Aquí se van a copiar todas las secuencias crudas)
Ir a la terminal de windows para copiar las secuencias crudas a la carpeta recién creada en el servidor
En la terminal local colocar la siguiente instrucción: scp -r -P 1967 D:/BEGIOAM/Secuencias_crudas/* celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Secuencias_crudas
verificar la cantidad de secuencias creadas con ls ls | wc -l (tienen que haber 20)---OK

una vez verificado se pasa al procesamiento

1- Visualizar un archivo zcat MED_ADR-Med139-S72.F.fq.gz | head -n 4
luego crear una carpeta para los fasqc: mkdir FASTQC_01
luego aplicar a cada secuencia el fastqc
fastqc *.fq.gz

Como resultado se van a generar archivos .html y .zip. (los voy a mover a la carpeta FASTQC_01)
Para mover a esta carpeta usar el comando mv *_fastqc.html *_fastqc.zip ../FASTQC_01/--OK

Crear la carpeta Trimmomatic_01 para copiar la secuencias analizadas y recortadas ahí

# Esto asume que tus archivos terminan en .F.fq.gz y .R.fq.gz
for infile in *.F.fq.gz; do

    # Extrae el nombre base quitando el .F.fq.gz
    base=$(basename "$infile" .F.fq.gz)
    
    echo "Procesando muestra: $base"
    
    trimmomatic PE -threads 10 -phred33 \
      "${base}.F.fq.gz" "${base}.R.fq.gz" \
      "../Trimmomatic_01/${base}.F_trim.fq.gz" "../Trimmomatic_01/${base}.F_untrim.fq.gz" \
      "../Trimmomatic_01/${base}.R_trim.fq.gz" "../Trimmomatic_01/${base}.R_untrim.fq.gz" \
      HEADCROP:12 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 subir a 100
done

Aplicar el FASTQC a estos archivos para ver como quedaron los recortes
crear una nueva carpeta Limpieza_secuencias_01 
ahora aplicar el FASTQC y mover los archivos a esta carpeta y luego al directorio para verlos
fastqc *trim.fq.gz  (en la carpeta de Trimmomatic_01)


Entrar a la carpeta Trimmomatic_01 y mover los .html y los .zip

mv *.html /home/celiabc/BEGIOAM/Secuencias_limpieza_01/
mv *.zip /home/celiabc/BEGIOAM/Secuencias_limpieza_01/

Mover los archivos.html del servidor a la computadora local para poder visualizarlos

scp -r -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Secuencias_limpieza_01/*.html D:/BEGIOAM/Secuencias_limpieza_01

Abrir los .html y verificar si mejoró la calidad. En caso de no mejorar, repetir a partir del paso 4, definiendo nuevos parámetros de limpieza.

Una vez obtenidos los .fastq con la calidad deseada, copiarlos del servidor a la computadora local para respaldarlos. Estos serán los datos que se emplearán en los análisis posteriores.

volver a correr trimmomatic incluyendo un filtro de calidad
revisar las nuevas secuencias
si todo está ok correrlas en FastP para eliminar las colas
hacer en ensamblaje
hacer el snp calling y hacer los filtros de los snps



######Aplicación nuevamente de TRimmomatic
#!/bin/bash
set -euo pipefail

OUTDIR="../Trimmomatic_02"
ADAPTERS="/ruta/a/adapters/TruSeq3-PE-2.fa"   # <-- ajusta a tu archivo de adaptadores

mkdir -p "$OUTDIR"

for infile in *.F.fq.gz; do

    # Extrae el nombre base quitando el .F.fq.gz
    base=$(basename "$infile" .F.fq.gz)

    # Verifica que exista el par reverso
    if [[ ! -f "${base}.R.fq.gz" ]]; then
        echo "AVISO: no se encontró ${base}.R.fq.gz, se omite la muestra." >&2
        continue
    fi

    echo "Procesando muestra: $base"

    trimmomatic PE -threads 10 -phred33 \
      "${base}.F.fq.gz" "${base}.R.fq.gz" \
      "${OUTDIR}/${base}.F_trim.fq.gz" "${OUTDIR}/${base}.F_untrim.fq.gz" \
      "${OUTDIR}/${base}.R_trim.fq.gz" "${OUTDIR}/${base}.R_untrim.fq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:30:10 \
      HEADCROP:12 \
      LEADING:20 \
      TRAILING:20 \
      SLIDINGWINDOW:4:20 \
      AVGQUAL:20 \
      MINLEN:100 \
      2> "${OUTDIR}/${base}.trimmomatic.log"

done

Aplicar FASTQC

fastqc *trim.fq.gz  (en la carpeta de Trimmomatic_02)


mv *.html /home/celiabc/BEGIOAM/Secuencias_limpieza_02/
mv *.zip /home/celiabc/BEGIOAM/Secuencias_limpieza_02/

Mover los archivos.html del servidor a la computadora local para poder visualizarlos

scp -r -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Secuencias_limpieza_02/*.html D:/BEGIOAM/Secuencias_limpieza_02

Abrir los .html y verificar si mejoró la calidad. En caso de no mejorar, repetir a partir del paso 4, definiendo nuevos parámetros de limpieza.

Una vez obtenidos los .fastq con la calidad deseada, copiarlos del servidor a la computadora local para respaldarlos. Estos serán los datos que se emplearán en los análisis posteriores.

Hacer el filtrado en fastp
instalarlo en el servidor
 # download the latest build
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

#!/bin/bash
set -euo pipefail

FASTP="$HOME/fastp"                  # ejecutable descargado con wget (como en tu imagen)
INDIR="../Trimmomatic_01"            # <-- carpeta con tus lecturas de entrada
OUTDIR="../fastp_polyGA"

mkdir -p "$OUTDIR"

for f in "${INDIR}"/*.F_trim.fq.gz; do

    base=$(basename "$f" .F_trim.fq.gz)
    r="${INDIR}/${base}.R_trim.fq.gz"

    if [[ ! -f "$r" ]]; then
        echo "AVISO: no se encontró $r, se omite." >&2
        continue
    fi

    echo "Procesando muestra: $base"

    "$FASTP" \
      -i "$f" -I "$r" \
      -o "${OUTDIR}/${base}.F_polyGA.fq.gz" \
      -O "${OUTDIR}/${base}.R_polyGA.fq.gz" \
      -A -Q \
      --trim_poly_g --poly_g_min_len 10 \
      --trim_poly_x --poly_x_min_len 10 \
      --length_required 50 \
      --thread 10 \
      --html "${OUTDIR}/${base}.fastp.html" \
      --json "${OUTDIR}/${base}.fastp.json" \
      2> "${OUTDIR}/${base}.fastp.log"

done

aplicar nuevamente trimmomatic

#!/bin/bash
set -euo pipefail

OUTDIR="../Trimmomatic_03"
ADAPTERS="/ruta/a/adapters/TruSeq3-PE-2.fa"   # <-- ajusta
CROPLEN=137                                    # <-- ajusta según la longitud real

mkdir -p "$OUTDIR"

for infile in *.F_polyGA.fq.gz; do

    base=$(basename "$infile" .F_polyGA.fq.gz)
    r="${base}.R_polyGA.fq.gz"

    if [[ ! -f "$r" ]]; then
        echo "AVISO: no se encontró $r, se omite la muestra." >&2
        continue
    fi

    echo "Procesando muestra: $base"

    trimmomatic PE -threads 10 -phred33 \
      "$infile" "$r" \
      "${OUTDIR}/${base}.F_trim.fq.gz" "${OUTDIR}/${base}.F_unpaired.fq.gz" \
      "${OUTDIR}/${base}.R_trim.fq.gz" "${OUTDIR}/${base}.R_unpaired.fq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:30:10 \
      HEADCROP:8 \
      CROP:${CROPLEN} \
      LEADING:20 \
      TRAILING:20 \
      SLIDINGWINDOW:4:20 \
      AVGQUAL:20 \
      MINLEN:100 \
      2> "${OUTDIR}/${base}.trimmomatic.log"

done

fastqc *trim.fq.gz  (en la carpeta de Trimmomatic_03)

mv *.html /home/celiabc/BEGIOAM/Secuencias_limpieza_03/
mv *.zip /home/celiabc/BEGIOAM/Secuencias_limpieza_03/

Mover los archivos.html del servidor a la computadora local para poder visualizarlos

scp -r -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Secuencias_limpieza_03/*.html D:/BEGIOAM/Secuencias_limpieza_03

#!/bin/bash
set -euo pipefail

INDIR="$HOME/BEGIOAM/Trimmomatic_02"        # <-- carpeta con los *_trim.fq.gz
OUTDIR="$HOME/BEGIOAM/Trimmomatic_04"
mkdir -p "$OUTDIR"

for f in "${INDIR}"/*.F_trim.fq.gz; do

    base=$(basename "$f" .F_trim.fq.gz)
    r="${INDIR}/${base}.R_trim.fq.gz"

    if [[ ! -f "$r" ]]; then
        echo "AVISO: no se encontró $r, se omite la muestra." >&2
        continue
    fi

    echo "Procesando muestra: $base"

    trimmomatic PE -threads 10 -phred33 \
      "$f" "$r" \
      "${OUTDIR}/${base}.F_trim.fq.gz" "${OUTDIR}/${base}.F_unpaired.fq.gz" \
      "${OUTDIR}/${base}.R_trim.fq.gz" "${OUTDIR}/${base}.R_unpaired.fq.gz" \
      CROP:125 \
      2> "${OUTDIR}/${base}.trimmomatic.log"

done

mv *.html /home/celiabc/BEGIOAM/Secuencias_limipieza_04/
mv *.zip /home/celiabc/BEGIOAM/Secuencias_limipieza_04/


####Ensamblar contra el genoma de referencia
copiar el genoma de referencia 

scp -r -P 1967 D:/BEGIOAM/genomev3/* celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Reference
verificar la cantidad de secuencias creadas con ls ls | wc -l 
-Indexar el genoma de referencia
# Crear un ambiente con todo lo necesario
conda create -n align -c conda-forge -c bioconda bwa-mem2 samtools fastp fastqc angsd bcftools
conda activate align

# Verificar
bwa-mem2 version

1- Se indexa

cd /home/celiabc/BEGIOAM/Reference
bwa index assembly3.fasta
ls assembly3.fasta*

2-Se alinea

cd /home/celiabc/BEGIOAM/Secuencias_limipieza_04
conda activate align
nohup bash alinear.sh > alinear.log 2>&1 &
sleep 60; tail -n 5 alinear.log

3-Se revisa la calidad del mapeo

conda activate align
cd /home/celiabc/BEGIOAM/Mapeo

for B in *.dedup.bam; do
  echo "== $B"
  samtools flagstat $B | grep -E "mapped \(|properly paired"
  samtools coverage $B | awk 'NR>1{l+=$3; c+=$7*$3} END{printf "Cobertura media: %.2fx\n", c/l}'
done

Busca >90% de lecturas mapeadas y una cobertura similar entre muestras. Si alguna muestra es muy baja (por ejemplo <1x o <70% mapeado), considera excluirla.



4-Buscar la lista de BAM
ls $PWD/*.dedup.bam > bam.list
wc -l bam.list      # debe dar 10

3. ANGSD: genotype likelihoods y SNPs

Con 10 muestras, -minInd 8 (80%). Ajusta -setMinDepth y -setMaxDepth según la cobertura que obtuviste (por ejemplo, mínimo 10 y máximo ~2x la profundidad total media).

angsd -bam bam.list -ref /home/celiabc/BEGIOAM/Reference/assembly3.fasta \
  -out chippurus \
  -GL 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
  -doCounts 1 -doPost 1 -doGeno 3 -postCutoff 0.95 \
  -minQ 20 -minMapQ 30 \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
  -minInd 8 -doBcf 1 -P 20



scp -r -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Mapeo/.tar.gz * D:/BEGIOAM/R




##análisis después de la poda

PCA y NAdmixture
cd /home/celiabc/BEGIOAM/Mapeo
mkdir -p PCA_podado NGSadmix_podado

2. Generar pruned.beagle.gz (si aún no lo hiciste)
bash
sed 's/:/_/' thin1kb.unlinked.id > unlinked.markers
awk 'NR==FNR{k[$1]=1; next} FNR==1 || ($1 in k)' unlinked.markers \
  <(zcat thin1kb.beagle.gz) | gzip > pruned.beagle.gz
zcat pruned.beagle.gz | wc -l     # debe dar 129049

PCA
conda activate popgen
cd /home/celiabc/BEGIOAM/Mapeo
pcangsd -b pruned.beagle.gz -o PCA_podado/pruned_pca -t 8 > PCA_podado/pcangsd.log 2>&1
ls -lh PCA_podado/


NAdmix

compactando archivos para pasarlos a la PC

cd /home/celiabc/BEGIOAM/Mapeo
tar czf podado_R.tar.gz PCA_podado/*.cov NGSadmix_podado/*.qopt NGSadmix_podado/*.log
tar tzf podado_R.tar.gz | wc -l

#copiar a la PC
scp -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Mapeo/podado_R.tar.gz D:/BEGIOAM/R/

















####6. SNPs calling en ipyrad

Crear la carpeta 04.ipyrad en /home/usuario/ Entrar a la carpeta y activar el ambiente conda para correr el programa.

conda activate ipyrad

Crear el archivo de parámetros que emplea ipyrad. Colocar un nombre informativo que permita identificar la corrida. En este ejemplo U90M90

ipyrad -n U90M90

Como resultado se obtendrá el archivo params-U90M90.txt. Este debe ser editado con el comando nano. Los parámetros deben ser modificados según el criterio del investigador. Se requieren varias corridas para definir el set de parámetros que optimiza el ensamblaje, permitiendo obtener el mayor número de SNPs de alta calidad, con bajo porcentaje de datos faltantes. Para más información consulta el manual de ipyrad

A continuación, algunos parámetros importantes a considerar:

    [1] [project_dir]: Especificar directorio de los archivos de salida.
    [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files. /ruta/*trim.fastq.gz
    [5] [assembly_method]: Assembly method (denovo, reference, etc)
    [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc. pair3rad
    [8] [restriction_overhang] colocar el sitio de corte de mis enzimas: AATTC, CTAGC, CTAGA (para EcoR1, Xba y NheI)
    [14] [clust_threshold]: Clustering threshold for de novo assembly. Colocar el umbral de similitud para la creación de clusters (en este ejemplo 0.90)
    [21] [min_samples_locus]: Min # samples per locus for output. Depende del tamaño de muestra y cuanto missing data acaptaré (ej. si tengo 60 muestras y solo aceptaré el 90% de missing data, en esta casilla debo colocar 54)
    [22] [max_SNPs_locus]: Número máximo de SNPs por locus. 0.1 = 10%
    [27] [output_formats]: Output formats (see docs). Colocar * para que genere todos los formatos posibles.

Una vez fijados los parámetros de la corrida modificando el params-U90M90, ya se puede correr ipyrad

nohup ipyrad -p params-U90M90.txt -s 1234567 -c 24 -f --MPI > output.txt &

Dependiendo de la catidad de datos, la corrida puede ser lenta. Recomiendo utilizar nohup para ir monitoreándola.

Una vez culminada la corrida, verificar los resultados guardados en la carpeta U90M90_outfiles, archivo U90M90_stats.txt. El archivo U90M90.vcf va a contener la matriz de SNPs.
7. Filtrado adicional con VCFtools

Crear la carpeta VCFtools en /home/usuario

Las opciones de filtrado en VCFtools pueden ser visualizadas con el comando man vcftools

Dependiendo del objetivo del estudio, se pueden hacer varias corridas. En este ejemplo, haremos análisis poblacionales clásicos (estructura y diversidad genética) y análisis de SNPs candidatos a selección. Por lo tanto, se construirán 2 sets de datos.

set de datos 1, análsis de selección

Objetivo: eliminar alelos con frecuencia menor al 4%, retener loci bialelélicos con profundidad mayor a 10 y máximo 10% de missing data.

vcftools --vcf U90M90.vcf --maf 0.04 --min-meanDP 10 --min-alleles 2 --max-alleles 2 --max-missing 0.9 --recode --out U90M90_maf0.04_DP10_MD0.9

set de datos 2, análisis de estructura y diversidad

Objetivo: A la salida anterior, aplicarle el filtro de desequilibrio de ligamiento.

vcftools --vcf U90M90_maf0.04_DP10_MD0.9.vcf --thin 1000 --recode --out U90M90_maf0.04_DP10_DL.vcf

Para ambos sets de datos, visualizar el archivo.log para ver cuantos sitios se mantienen luego del filtrado.

En este punto, ya están listos los sets de datos a emplear en análisis de estructura, diversidad y selección. Revisar el nombre de las muestras en los vcf resultantes y modificarlos si es requerido.
Descargar el vcf

En la terminal de mi compu scp -r -P 1967 massiel.alfonsoglez@132.248.15.41:/home/massiel.alfonsoglez/curso/04.ipyrad/U90M90_outfiles/U90M90_maf0.04_DP10_MD0.9.vcf /home/alex-dll/Escritorio



 
