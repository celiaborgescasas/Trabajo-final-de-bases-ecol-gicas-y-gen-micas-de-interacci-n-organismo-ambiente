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

 
