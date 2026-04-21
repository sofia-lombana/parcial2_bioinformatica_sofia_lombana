# parcial2_bioinformatica_sofia_lombana

# **Parcial**

## **Punto 1**
**a):**
- Per base sequence quality - PASS - buena calidad
- Per sequence quality scores - PASS - lecturas confiables
- Adapter content - PASS - no hay adaptadores
- Overrepresented sequences - PASS - no contaminación evidente
- Per base sequence content - FAIL - hay diferencias en proporciones A, T, G, C al inicio

[FastQC 1](https://github.com/sofia-lombana/parcial2_bioinformatica_sofia_lombana/blob/main/fastqc_1.jpg)

La calidad de las secuencias es alta en la mayor parte del read (Q>30), lo que indica datos confiables. Se observa una leve disminución hacia el final, algo normal en secuencias Illumina, pero sin afectar significativamente el análisis.

[FastQC 2](https://github.com/sofia-lombana/parcial2_bioinformatica_sofia_lombana/blob/main/fastqc_2.jpg)

La calidad de las secuencias es alta en la mayor parte del read (Q>30), aunque se observa una ligera disminución hacia el final. Esta caída es esperada en datos Illumina y no compromete significativamente la calidad global de las secuencias.

Codigo:
1. ```nano fastqc.sh```

2. ```
#!/bin/bash
#SBATCH -p normal #Particion (cola)
#SBATCH -N 1 # Numero de nodos
#SBATCH -n 8 # Numero de nucleos
#SBATCH -t 0-00:20 # limite de tiempo (D-HH:MM)
#SBATCH -o fastqc.out # Salida STDOUT
#SBATCH -e fastqc.err # Salida STDERR
#SBATCH --mail-user=sofia.lombana@urosario.edu.co #Direccion e-mail a donde notificar el estado del trabajo
#SBATCH --mail-type=ALL #Especifica que eventos notificar al correo (ALL, BEGIN, END, REQUEUE, FAIL)
module load fastqc
fastqc SRR8528336_1.fastq.gz SRR8528336_2.fastq.gz```

3. ```sbatch fastqc.sh```

**b):**

El ensamblaje con SPAdes se realizó correctamente, ya que se generaron los archivos contigs.fasta y scaffolds.fasta sin errores, indicando que el proceso fue exitoso.

Codigo:
1.  ```nano spades.sh```
2. ```
#!/bin/bash
#SBATCH -p normal #Particion (cola)
#SBATCH -N 1 # Numero de nodos
#SBATCH -n 8 # Numero de nucleos
#SBATCH -t 0-00:20 # limite de tiempo (D-HH:MM)
#SBATCH -o spades.out # Salida STDOUT
#SBATCH -e spades.err # Salida STDERR
#SBATCH --mail-user=sofia.lombana@urosario.edu.co #Direccion e-mail a donde notificar el estado del trabajo
#SBATCH --mail-type=ALL #Especifica que eventos notificar al correo (ALL, BEGIN, END, REQUEUE, FAIL)
module load spades
spades.py --careful -1 SRR8528336_1.fastq.gz -2 SRR8528336_2.fastq.gz -o spades_ensamblaje
```
3. ```sbatch spades.sh```

**c):**

Métrica | Valor
:----------:|:--------------:
Tamaño total ensamblado    | 2,100,634 bp
N50      | 1,178 bp
L50     | 507

El ensamblaje es altamente fragmentado y de baja calidad, con un N50 bajo y un L50 alto, lo que indica que el genoma está dividido en muchos contigs pequeños y no se logró una reconstrucción continua.

Codigo:
1. ```/datacnmat01/ciencias/appsbio/conda/pkgs/quast-5.2.0-py39pl5321h4e691d4_3/opt/quast-5.2.0/quast.py scaffolds.fasta -o quast_results```
2. ```QUAST
Assembly                    scaffolds
# contigs (>= 0 bp)         5426
# contigs (>= 1000 bp)      632
# contigs (>= 5000 bp)      7
# contigs (>= 10000 bp)     1
# contigs (>= 25000 bp)     0
# contigs (>= 50000 bp)     0
Total length (>= 0 bp)      3089161
Total length (>= 1000 bp)   1185860
Total length (>= 5000 bp)   56955
Total length (>= 10000 bp)  15601
Total length (>= 25000 bp)  0
Total length (>= 50000 bp)  0
# contigs                   1970
Largest contig              15601
Total length                2100634
GC (%)                      41.76
N50                         1178
N90                         573
auN                         1701.7
L50                         507
L90                         1578
# N's per 100 kbp           102.68 
```

## **Punto 2**
**a):**

Codigo:
1. ```makeblastdb -in contigs.fasta -dbtype nucl -out contigs_db```
2. ```blastn -query sec.fasta -db contigs_db -out blast_obp.txt -outfmt 7```

**b):**

Codigo:
1. ```sed -E 's/^>(\w+\.\w+) (\w+) (\w+).*/>\1_\2_\3/' sec.fasta > sec_final.fasta```
2. ```sed -i -E 's/^>(\w+\.\w+) (\w+:) (\w+) (\w+).*/>\1_\3_\4/' sec_final.fasta```

**c):**

Codigo:
```
blastn -query sec_final.fasta \
-db contigs_db \
-outfmt 7 \
-out blast_obp:final.txt \
-num_threads 1
```

Se realizó una búsqueda BLAST utilizando las secuencias del archivo `sec_final.fasta` como query contra la base de datos local construida a partir del ensamblaje (`contigs.fasta`). 

No se encontraron hits significativos en el archivo de salida (`blast_final.txt`), lo que indica que no fue posible identificar un ortólogo del gen ommochrome-binding protein (OBP) en el ensamblaje obtenido.

Este resultado puede explicarse por la alta fragmentación del ensamblaje, evidenciada previamente por un valor bajo de N50 y un valor alto de L50, lo que dificulta la reconstrucción de genes completos. Adicionalmente, la ausencia de coincidencias podría estar relacionada con diferencias evolutivas entre las secuencias de referencia y el organismo analizado.

En conjunto, estos resultados sugieren que la calidad y contigüidad del ensamblaje limitan la detección del gen de interés.

## **Punto 3**

**a):**

Codigo:
```
> alt <- read.csv("df_altitud.csv")
> head (alt)
  X Altitud    Grupo   Abundancia
1 1     500     Rojo 4.108510e-04
2 2     500 Amarillo 3.427287e-06
3 3     500   Blanco 9.503885e-15
4 4     510     Rojo 4.277928e-04
5 5     510 Amarillo 3.689888e-06
6 6     510   Blanco 1.165079e-14
> library(ggplot2)
> pdf("altitud_morfos.pdf", width = 8, height = 5)
> ggplot(alt, aes(x = Altitud, y = Abundancia, fill = Grupo)) +
+   geom_area(position = "fill") +
+   scale_fill_manual(values = c("Blanco" = "white",
+                                "Amarillo" = "yellow",
+                                "Rojo" = "red")) +
+   labs(title = "Distribución de morfos según altitud",
+        x = "Altitud",
+        y = "Abundancia relativa",
+        fill = "Fenotipo") +
+   theme_bw()
> dev.off ()
```
[Grafica](https://github.com/sofia-lombana/parcial2_bioinformatica_sofia_lombana/blob/main/altitud_morfos.pdf)


**b):**

Codigo:
```
> temp <- read.csv("df_temperatura.csv")
> head (temp)
  X Temperatura Supervivencia Grupo
1 1          22    0.04997452  Rojo
2 2          22    0.12762019  Rojo
3 3          29    0.66976309  Rojo
4 4          21    0.00000000  Rojo
5 5          25    0.38718651  Rojo
6 6          30    0.96059158  Rojo
> library(ggplot2)
> pdf("temperatura_supervivencia.pdf", width = 8, height = 5)
> ggplot(temp, aes(x = Temperatura, y = Supervivencia, color = Grupo)) +
+   geom_point(size = 2) +
+   scale_color_manual(values = c("Blanco" = "black",
+                                 "Amarillo" = "orange",
+                                 "Rojo" = "red")) +
+   labs(title = "Relación entre temperatura y probabilidad de supervivencia",
+        x = "Temperatura",
+        y = "Probabilidad de supervivencia",
+        color = "Fenotipo") +
+   theme_bw()
> dev.off()
```
[Grafica](https://github.com/sofia-lombana/parcial2_bioinformatica_sofia_lombana/blob/main/temperatura_supervivencia.pdf)
