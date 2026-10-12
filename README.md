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

####Ensamblar contra el genoma de referencia
esto se hace con ipayrad

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



 
