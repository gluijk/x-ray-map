# Procesado de mapas raster con R
# www.overfitting.net
# https://www.overfitting.net/2020/10/procesado-de-mapas-raster-con-r.html

library(data.table)  # fread()
library(tiff)


# DIBUJANDO MAPA 3D DE ELEVACIONES DESDE DATOS RASTER

# Centro de Descargas del Centro Nacional de Información Geográfica
# Modelos de elevaciones en formato raster MDT25 (resolución rejilla=25m)
# URL: http://centrodedescargas.cnig.es/CentroDescargas/index.jsp

# Leemos y procesamos datos raster
# 7 cuadrantes cubriendo la isla de Tenerife
# Cotas en m, resolución rejilla=25m
tenerife1=data.matrix(
    fread("PNOA_MDT25_REGCAN95_HU28_1088_LID.txt", sep=" ", dec="."))
tenerife2=data.matrix(
    fread("PNOA_MDT25_REGCAN95_HU28_1089_LID.txt", sep=" ", dec="."))
tenerife3=data.matrix(
    fread("PNOA_MDT25_REGCAN95_HU28_1091_LID.txt", sep=" ", dec="."))
tenerife4=data.matrix(
    fread("PNOA_MDT25_REGCAN95_HU28_1092_LID.txt", sep=" ", dec="."))
tenerife5=data.matrix(
    fread("PNOA_MDT25_REGCAN95_HU28_1096_LID.txt", sep=" ", dec="."))
tenerife6=data.matrix(
    fread("PNOA_MDT25_REGCAN95_HU28_1097_LID.txt", sep=" ", dec="."))
tenerife7=data.matrix(
    fread("PNOA_MDT25_REGCAN95_HU28_1102_LID.txt", sep=" ", dec="."))

# Crop sea areas to 0
tenerife1[tenerife1<0]=0
tenerife2[tenerife2<0]=0
tenerife3[tenerife3<0]=0
tenerife4[tenerife4<0]=0
tenerife5[tenerife5<0]=0
tenerife6[tenerife6<0]=0
tenerife7[tenerife7<0]=0
MAXIMO=max(tenerife3)  # Teide: 3710.062m

writeTIFF(tenerife1/MAXIMO, "tenerife1.tif", bits.per.sample=16, compression="LZW")
writeTIFF(tenerife2/MAXIMO, "tenerife2.tif", bits.per.sample=16, compression="LZW")
writeTIFF(tenerife3/MAXIMO, "tenerife3.tif", bits.per.sample=16, compression="LZW")
writeTIFF(tenerife4/MAXIMO, "tenerife4.tif", bits.per.sample=16, compression="LZW")
writeTIFF(tenerife5/MAXIMO, "tenerife5.tif", bits.per.sample=16, compression="LZW")
writeTIFF(tenerife6/MAXIMO, "tenerife6.tif", bits.per.sample=16, compression="LZW")
writeTIFF(tenerife7/MAXIMO, "tenerife7.tif", bits.per.sample=16, compression="LZW")

# Put all 7 images together manually in Photoshop -> "tenerifecomposite.tif"
tenerife=readTIFF("tenerifecomposite.tif")*MAXIMO

hist(tenerife[tenerife>0], breaks=800)
