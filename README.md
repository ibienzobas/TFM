# TFM
"Identificación de alelos de genes de síntesis de volátiles en el genoma de melón".
Conjunto de resultados generados para el TFM, asi como la herramienta bioinformática que se ha creado.

Estructura de resultados: 

Para cada enzima involucrada en la síntesis de VOCs estudiada: 

- Archivo ".genes": con los posibles genes en el genoma de melón que codifiquen dicha enzima

- 8 archivos ".mutations": uno por cada variedad de melón estudiada. Contienen las mutaciones encontradas en los posibles genes codificantes de la enzima (archivo ".genes"), para cada variedad de melón.

- Archivo "unique_mutations": con las mutaciones únicas encontradas en la variedad de melón seleccionada (no aparecen en el resto de variedades de melón), en los posibles genes de sínteis de la enzima. Este archivo solo se encuentra en para las 4 enzimas y estudiadas mas en detalle en el TFM (y las variedades correspondientes): (S)-1-Feniletanol deshidrogenasa (para el compuesto acetofenona), ZCD (ciclocitral), CCD1 (ionona), nonanal (LOX1). Pero puede ser generado para la enzima y variedad de melón que se desee, siempre que el script se lance con el argumento "-m" (modo mutaciones únicas).











Opciones de lanzamiento del script: 

1- Modo generación de datos: Genera los archivos ".genes" y ".mutations" al suministrarle como input el directorio de trabajo (donde genera dichos archivos), el directorio donde se encuentran las variedades de melón estudiadas (VCFs), el genoma anotado de melón y los archivos de salida del tBLASTn. 

Formato de lanzamiento: 

'script.py  -b < tBLASTn directory> -v <variants directory> -g <genome file> -o <output directory(working directory)>'




2- Modo mutaciones únicas: Genera el archivo "unique_mutations" para la variedad de melón y la enzima seleccionadas. Necesita del directorio donde se han generado los resultados (directorio de trabajo del modo generación de datos). El propio script lee las enzimas y variedades disponibles y hace elegir una de cada. Debe añadir el argumento "-m" para activar este modo.

Formato de lanzamiento:

'script.py -m -d <results directory>'



Para consultar los formatos de lanzamiento del script siempre puede lanzarlo junto con el argumento "-h".