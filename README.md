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

angsd -bam bam.list -ref /home/celiabc/BEGIOAM/Reference/assembly3.fasta \
  -out chippurus \
  -GL 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
  -doCounts 1 -doDepth 1 -dumpCounts 2 \
  -setMinDepth 10 -setMaxDepth 50 \
  -minQ 20 -minMapQ 30 \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
  -minInd 8 \
  -doGlf 2 -doBcf 1 -doPost 1 -P 8



scp -r -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Mapeo/.tar.gz * D:/BEGIOAM/R

#POda
Paso 1. Adelgazar: un SNP cada ≥ 1 kb
bash
cd /home/celiabc/BEGIOAM/Mapeo
zcat chippurus.mafs.gz | tail -n +2 | cut -f1,2 | \
  awk '{ if ($1!=c || $2-last>=1000) {print; c=$1; last=$2} }' > thin1kb.pos
wc -l thin1kb.pos        # 379,547

Paso 2. Extraer esos sitios del archivo de likelihoods
bash
awk 'NR==FNR{k[$1"_"$2]=1; next} FNR==1 || ($1 in k)' thin1kb.pos \
  <(zcat chippurus.beagle.gz) | gzip > thin1kb.beagle.gz
zcat thin1kb.beagle.gz | wc -l      # 379,548 (sitios + encabezado)

Paso 3. Calcular el LD con ngsLD
bash
conda activate popgen
NS=$(wc -l < thin1kb.pos)
nohup ngsLD --geno thin1kb.beagle.gz --pos thin1kb.pos --probs \
  --n_ind 10 --n_sites $NS --max_kb_dist 10 --min_maf 0.05 \
  --n_threads 4 --out thin1kb.ld > ngsld.log 2>&1 &

Paso 4. Podar con prune_graph
bash
export PATH=$HOME/.cargo/bin:$PATH
prune_graph --in thin1kb.ld --header \
  --weight r2 \
  --filter "dist <= 10000 && r2 >= 0.2" \
  --keep-heavy \
  --out thin1kb.unlinked.id \
  -n 4
wc -l thin1kb.unlinked.id      # 129,048

Paso 5. Crear el archivo de likelihoods podado
bash
sed 's/:/_/' thin1kb.unlinked.id > unlinked.markers
awk 'NR==FNR{k[$1]=1; next} FNR==1 || ($1 in k)' unlinked.markers \
  <(zcat thin1kb.beagle.gz) | gzip > pruned.beagle.gz
zcat pruned.beagle.gz | wc -l        # 129,049

Resumen de números
Etapa	SNPs
ANGSD (todos)	3,504,603
Tras adelgazar (≥ 1 kb)	379,547
Tras prune_graph (r² < 0.2, ≤ 10 kb)	129,048

###Siguientes análisis

PCA
#sin control de ligamiento
conda activate popgen
cd /home/celiabc/BEGIOAM/Mapeo
pcangsd -b chippurus.beagle.gz -o chippurus_pca -t 8
mkdir -p PCA_podado
pcangsd -b pruned.beagle.gz -o PCA_podado/pruned_pca -t 8 > PCA_podado/pcangsd.log 2>&1

#con control de ligamiento
conda activate popgen
cd /home/celiabc/BEGIOAM/Mapeo
pcangsd -b pruned.beagle.gz -o PCA_podado/pruned_pca -t 8 > PCA_podado/pcangsd.log 2>&1
ls -lh PCA_podado/

###Visualización en R
library(ggplot2); library(ggrepel)
setwd("D:/BEGIOAM/R")

# Orden de bam.list (el mismo del .cov)
inds <- c("MED_ADR-Med139-S72","MED_ADR-Med146-S35","MED_ADR-Med146-S88",
          "MED_CHI-CH18C-S53","MED_CHI-CH9C-S41",
          "ORA_DAK-DAK3-S46","ORA_DAK-DAK3-S50","ORA_DAK-DAK4-S62",
          "ORA_DAK-DAK5-S74","ORA_DAK-DAK6-S86")
pob <- factor(sub("-.*", "", inds),
              levels = c("MED_ADR", "MED_CHI", "ORA_DAK"),
              labels = c("Adriático", "Chipre", "Dakar"))
cols <- c("Adriático" = "#4477AA", "Chipre" = "#EE6677", "Dakar" = "#228833")

# Matriz de covarianza -> componentes principales
cov_mat <- as.matrix(read.table("pruned_pca.cov"))   # sin poda: "chippurus_pca.cov"
e   <- eigen(cov_mat)
var <- round(100 * e$values / sum(e$values), 1)

pca <- data.frame(ind = sub("^[A-Z_]+-", "", inds), pob = pob,
                  PC1 = e$vectors[, 1], PC2 = e$vectors[, 2], PC3 = e$vectors[, 3])

p <- ggplot(pca, aes(PC1, PC2, color = pob, label = ind)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = cols) +
  labs(title = "PCA con poda por LD (129,048 SNPs)",
       x = paste0("PC1 (", var[1], "%)"),
       y = paste0("PC2 (", var[2], "%)"), color = "Población") +
  theme_bw()

print(p)
ggsave("PCA_con_poda.png", p, width = 7, height = 5, dpi = 300)
NAdmix

compactando archivos para pasarlos a la PC

cd /home/celiabc/BEGIOAM/Mapeo
tar czf podado_R.tar.gz PCA_podado/*.cov NGSadmix_podado/*.qopt NGSadmix_podado/*.log
tar tzf podado_R.tar.gz | wc -l

#copiar a la PC
scp -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Mapeo/podado_R.tar.gz D:/BEGIOAM/R/


revisar gráficos en R

scp -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Mapeo/metadata.tsv D:/BEGIOAM/R/


###Siguiente análisis: ancestría con NGSadmix (SNPs podados)

Después del PCA corriste NGSadmix sobre pruned.beagle.gz, y luego lo graficaste en R.

Terminal (servidor, ambiente align)

conda activate align
cd /home/celiabc/BEGIOAM/Mapeo
mkdir -p NGSadmix_podado

nohup bash -c '
for K in 1 2 3 4; do
  for R in 1 2 3 4 5; do
    NGSadmix -likes pruned.beagle.gz -K $K -P 8 -seed $R -o NGSadmix_podado/pruned_K${K}_r${R} > /dev/null
  done
done
echo TERMINADO
' > NGSadmix_podado/ngsadmix.log 2>&1 &

-likes: entrada de likelihoods (beagle).
-K 1..4: número de grupos ancestrales probados.
-seed 1..5: 5 réplicas por K, porque el algoritmo puede caer en óptimos locales.
Salidas: .qopt (ancestría por individuo), .log (verosimilitud), .fopt.gz.


#Elegir K (verosimilitud) y ancestría

grep -H "best like" NGSadmix_podado/pruned_K*_r*.log | sed 's#NGSadmix_podado/##' | sort -V
paste metadata.tsv NGSadmix_podado/pruned_K2_r1.qopt | column -t

#Empaquetar y copiar a tu PC:
tar czf podado_R.tar.gz PCA_podado/*.cov NGSadmix_podado/*.qopt NGSadmix_podado/*.log
# en la terminal de tu PC
scp -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Mapeo/podado_R.tar.gz D:/BEGIOAM/R/


#En R
setwd("D:/BEGIOAM/R")
library(ggplot2); library(tidyr); library(dplyr)

inds <- c("MED_ADR-Med139-S72","MED_ADR-Med146-S35","MED_ADR-Med146-S88",
          "MED_CHI-CH18C-S53","MED_CHI-CH9C-S41",
          "ORA_DAK-DAK3-S46","ORA_DAK-DAK3-S50","ORA_DAK-DAK4-S62",
          "ORA_DAK-DAK5-S74","ORA_DAK-DAK6-S86")
meta <- data.frame(muestra = sub("^[A-Z_]+-", "", inds), pob = sub("-.*", "", inds))

# Verosimilitud por corrida
logs <- list.files("NGSadmix_podado", pattern = "^pruned_K[0-9]+_r[0-9]+\\.log$", full.names = TRUE)
get_ll <- function(f) {
  l <- grep("best like", readLines(f), value = TRUE)
  as.numeric(sub(".*best like=([-0-9.eE+]+).*", "\\1", l[length(l)]))
}
res <- data.frame(archivo = basename(logs),
                  K = as.numeric(sub(".*_K([0-9]+)_r.*", "\\1", logs)),
                  rep = as.numeric(sub(".*_r([0-9]+)\\.log", "\\1", logs)),
                  loglike = sapply(logs, get_ll))
print(res[order(res$K, res$rep), ], row.names = FALSE)
print(res %>% group_by(K) %>% summarise(media = mean(loglike), maximo = max(loglike), sd = sd(loglike)))

ggplot(res, aes(factor(K), loglike)) + geom_point(size = 3, alpha = .7) +
  labs(x = "K", y = "Log-verosimilitud") + theme_bw()

# Barras de ancestría (mejor réplica de cada K)
colores <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
graficar_K <- function(k) {
  mejor <- res %>% filter(K == k) %>% slice_max(loglike, n = 1, with_ties = FALSE)
  q <- read.table(file.path("NGSadmix_podado", sub("\\.log$", ".qopt", mejor$archivo)))
  colnames(q) <- paste0("Grupo", 1:k)
  d <- pivot_longer(cbind(meta, q), starts_with("Grupo"), names_to = "grupo", values_to = "prop")
  ggplot(d, aes(muestra, prop, fill = grupo)) + geom_col(width = 1) +
    facet_grid(~pob, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = colores) +
    labs(x = NULL, y = "Ancestría", title = paste0("K = ", k)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
}
graficar_K(2)


#Terminal (servidor, ambiente align)

Todo se lanzó con un script, FST_podado.sh:
cd /home/celiabc/BEGIOAM/Mapeo
cat > FST_podado.sh << 'EOF'
#!/bin/bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate align
cd /home/celiabc/BEGIOAM/Mapeo
REF=/home/celiabc/BEGIOAM/Reference/assembly3.fasta
FILT="-minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1"
mkdir -p FST_podado

# 1) Sitios podados en formato ANGSD (contig<TAB>pos) e índice
sed 's/:/\t/' thin1kb.unlinked.id | sort -k1,1 -k2,2n > FST_podado/pruned.sites
angsd sites index FST_podado/pruned.sites
[ -f MED.list ] || cat MED_ADR.list MED_CHI.list > MED.list

# 2) SAF por población, solo en los sitios podados
for P in MED ORA_DAK MED_ADR MED_CHI; do
  N=$(wc -l < $P.list)
  MI=$(( N > 2 ? N-1 : N ))
  angsd -bam $P.list -ref $REF -anc $REF -sites FST_podado/pruned.sites \
    -out FST_podado/saf_$P -GL 2 -doSaf 1 -doCounts 1 $FILT \
    -minInd $MI -P 4 > FST_podado/saf_$P.log 2>&1 &
done
wait

# 3) FST por pares
cd FST_podado
for PAIR in "MED ORA_DAK" "MED_ADR ORA_DAK" "MED_CHI ORA_DAK" "MED_ADR MED_CHI"; do
  set -- $PAIR
  realSFS saf_$1.saf.idx saf_$2.saf.idx -P 8 > $1_$2.sfs
  realSFS fst index saf_$1.saf.idx saf_$2.saf.idx -sfs $1_$2.sfs -fstout fst_$1_$2
  echo "== $1 vs $2"
  realSFS fst stats fst_$1_$2.fst.idx
  realSFS fst print fst_$1_$2.fst.idx > persite_$1_$2.txt
done | tee fst_global.txt
echo TERMINADO
EOF
chmod +x FST_podado.sh
nohup ./FST_podado.sh > FST_podado.log 2>&1 &

#Seguimiento de los resultados
tail -n 2 FST_podado.log          # TERMINADO
cat FST_podado/fst_global.txt     # FST no ponderado y ponderado

#Empaquetar para copiar a la PC
tar czf FST_podado_R.tar.gz FST_podado/fst_global.txt FST_podado/persite_*.txt
#Copiar a la PC
scp -P 1967 celiabc@132.248.15.41:/home/celiabc/BEGIOAM/Mapeo/FST_podado_R.tar.gz D:/BEGIOAM/R/


#En R
library(ggplot2)
setwd("D:/BEGIOAM/R/FST_podado")

pares <- c("MED_ORA_DAK", "MED_ADR_ORA_DAK", "MED_CHI_ORA_DAK", "MED_ADR_MED_CHI")

res <- lapply(pares, function(p) {
  d <- read.table(paste0("persite_", p, ".txt"), col.names = c("chr", "pos", "A", "B"))
  d <- d[d$B > 0, ]
  est <- sum(d$A) / sum(d$B)                 # FST ponderado = sum(A)/sum(B)
  set.seed(1)
  bs <- replicate(1000, { i <- sample(nrow(d), replace = TRUE); sum(d$A[i]) / sum(d$B[i]) })
  data.frame(par = p, nSNP = nrow(d), FST = round(est, 4),
             IC_inf = round(quantile(bs, .025), 4),
             IC_sup = round(quantile(bs, .975), 4))
})
tabla <- do.call(rbind, res); rownames(tabla) <- NULL
print(tabla)
write.csv(tabla, "tabla_FST_podado.csv", row.names = FALSE)

ggplot(tabla, aes(par, FST)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup), width = .15) +
  labs(x = NULL, y = "FST ponderado (IC 95%)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("FST_podado.png", width = 6, height = 4, dpi = 300)

