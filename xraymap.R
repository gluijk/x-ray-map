# X-Ray animation from a Tenerife DEM
# www.overfitting.net
# https://www.overfitting.net/2024/10/radiografia-de-tenerife-con-r.html

library(data.table)  # fread()
library(terra)  # build blur and resample functions
library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's


# Generic array resample function
# works both for matrix (grayscale images) or 3-channel arrays (colour images)
arrayresample=function(img, DIMX, DIMY, method='bilinear') {
    require(terra)
    
    raster=rast(img)
    rasterrs=rast(nrows=DIMY, ncols=DIMX, extent=ext(raster))
    rasterrs=resample(raster, rasterrs, method=method)
    
    if (is.matrix(img)) return (matrix(as.array(rasterrs), nrow=nrow(rasterrs)))
    else return (as.array(rasterrs))  # convert back to matrix/array
}


# Blur
# https://stackoverflow.com/questions/70429190/how-can-i-perform-neighborhood-analysis-in-terra-or-raster-and-keep-the-same-na
arrayblur=function(img, radius=11) {
    # radius: radius of the circular averaging window
    require(terra)
    
    # Build circular kernel
    D=2*radius+1  # D will always be an odd number as required by focal()
    w=matrix(1, nrow=D, ncol=D)
    w[(row(w)-(radius+1))^2 + (col(w)-(radius+1))^2 >= (radius+1)^2]=NA
    writePNG(w, "blurkernel.png")
    
    raster=rast(img)  # process as raster
    rasterblur=focal(raster, w=w, fun='mean', na.rm=TRUE, na.policy='omit')
    
    if (is.matrix(img)) return (matrix(as.array(rasterblur), nrow=nrow(rasterblur)))
    else return (as.array(rasterblur))  # convert back to matrix/array
}


# Hillshade calculation
hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}


# Shadow projection calculation
shadowmap=function(DEM, dx=25, dlight=c(0, 2, 3)) {
    # shadowsmap() inputs DEM data and outputs shadows projection
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    
    # Loop version: 2.22min (FASTER than the vectorized version)
    
    # Check for invalid input parameters
    if (dlight[3]<=0) {
        print(paste0("ERROR: Z light source must be positive (Z=", dlight[3], ")"))
        return (-1)
    } else if (dlight[1] & dlight[2]) {
        print(paste0("ERROR: X or Y light source must be 0 (X=", dlight[1], ", Y=", dlight[2],  ")"))
        return (-1)
    }
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    # Turn any light source into East light source setting dlightY
    if (!dlight[1] & !dlight[2]) {
        print(paste0("WARNING: zenithal Z light source, no shadows (X=", dlight[1], ", Y=", dlight[2], ")"))
        return (DEM*0+1)  # return white matrix
    } else if (dlight[2]<0) {  # West light source
        dlightY=-dlight[2]
        DEM=DEM[,ncol(DEM):1]  # transpose cols
    } else if (dlight[1]>0) {  # South light source
        dlightY=dlight[1]
        DEM=t(DEM)  # transpose all
    } else if (dlight[1]<0) {  # North light source
        dlightY=-dlight[1]
        DEM=t(DEM)  # transpose all
        DEM=DEM[,ncol(DEM):1]  # transpose cols
    } else {  # East light source (standard case)
        dlightY=dlight[2]
    }
    dlightZ=dlight[3]
    
    # DIMY and DIMX change if South/North light source
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    
    shadows=DEM*0  # no shadows at the beginning
    for (y in 1:(DIMX-1)) {  # col DIMX excluded since it cannot get shadow from anyone
        LEN=DIMX-y+1
        dh=(LEN-1)*dx*dlightZ/dlightY
        for (x in 1:DIMY) {
            LIGHTPATH=seq(from=DEM[x,y], to=DEM[x,y]+dh, length.out=LEN)
            DELTA=DEM[x, y:DIMX]-LIGHTPATH
            
            # 4 shadow styles:
            if (max(DELTA)>0) shadows[x,y]=1  # style='hard': shadow=1 if some elevation protrudes above light path
            # shadows[x,y]=shadows[x,y] + length(DELTA[DELTA>0])  # style='width': thickness of protrusion
            # shadows[x,y]=shadows[x,y] + max(0,DELTA[DELTA>0])  # style='max': highest protrusion
            # shadows[x,y]=shadows[x,y] + sum(DELTA[DELTA>0])  # style='mixed': width + height of protrusion
        }
    }
    
    # Undo transformations if applied
    if (dlight[2]<0) {  # West light source
        shadows=shadows[,ncol(shadows):1]  # transpose cols
    } else if (dlight[1]>0) {  # South light source
        shadows=t(shadows)  # transpose all
    } else if (dlight[1]<0) {  # North light source
        shadows=shadows[,ncol(shadows):1]  # transpose cols
        shadows=t(shadows)  # transpose all
    }
    
    return(1-shadows/max(shadows))  # normalize 0..1, shadow=black
    
    
    # Vectorized version: 3.36min (SLOWER than the loop version)
    
    # First we vectorize the needed linear interpolation
    # seq.vectorized=Vectorize(seq.default, vectorize.args=c("from", "to"))
    # 
    # shadows=DEM*0  # 0=no shadow, 1=shadow (no shadow at the beginning)
    # for (y in 1:(DIMX-1)) {  # col DIMX excluded since it cannot get shadow from anyone
    #     LEN=DIMX-y+1
    #     dh=(LEN-1)*dx*dlightZ/dlightY
    #     LIGHTPATH=t(seq.vectorized(from=DEM[,y], to=DEM[,y]+dh, length.out=LEN))
    #     DELTA=DEM[, y:DIMX]-LIGHTPATH
    #     shadows[,y]=apply(DELTA, MARGIN=1, FUN=max)  # max of each row
    # }
    # shadows[shadows<=0]=0
    # shadows[shadows>0]=1
}


# Contours calculation
contour=function(DEM, stroke=1) {
    # contour() calculates the contours of any colour change
    #
    # stroke: line width (pixels)
    
    DIMY=nrow(DEM)    
    DIMX=ncol(DEM)
    
    # Calculate outline map from solid map
    outline=DEM*0
    
    # 1 pixel stroke outline
    outline[2:(DIMY-1), 2:(DIMX-1)]=
        abs(DEM[1:(DIMY-2), 2:(DIMX-1)] -
                DEM[2:(DIMY-1), 2:(DIMX-1)]) +
        abs(DEM[2:(DIMY-1), 1:(DIMX-2)] -
                DEM[2:(DIMY-1), 2:(DIMX-1)])
    
    # 2 pixel stroke outline
    if (stroke==2 | stroke>3) {
        outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)] +
            outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
    }
    
    # 3 pixel stroke outline
    if (stroke>2) {
        outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)] +
            outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)] +
            outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]
    }
    
    outline[outline!=0]=1
    return(outline)
}


# Mask overlap calculation
overlap=function(DIMX, DIMY, POSX, POSY, dimx, dimy) {
    # Function generated using ChatGPT
    # prompt: "Write a R function that provides the top left corner
    # and bottom right corner of the overlapping between a rectangle
    # defined by coordinates (1,1) and (DIMX, DIMY) and another rectangle
    # of dimensions dimx by dimy, whose top left corner is positioned
    # at coordinates (POSX, POSY)"
    
    # Top-left corner of rectangle 1
    r1_x1=1
    r1_y1=1
    # Bottom-right corner of rectangle 1
    r1_x2=DIMX
    r1_y2=DIMY
    
    # Top-left corner of rectangle 2
    r2_x1=POSX
    r2_y1=POSY
    # Bottom-right corner of rectangle 2
    r2_x2=POSX + dimx - 1  # -1 added afterwards to work with
    r2_y2=POSY + dimy - 1  # integer matrix positions
    
    # Calculate top-left corner of the overlap
    overlap_x1=max(r1_x1, r2_x1)
    overlap_y1=max(r1_y1, r2_y1)
    
    # Calculate bottom-right corner of the overlap
    overlap_x2=min(r1_x2, r2_x2)
    overlap_y2=min(r1_y2, r2_y2)
    
    # Check if there is no overlap
    if (overlap_x1>overlap_x2 || overlap_y1>overlap_y2) {
        return(NULL)  # No overlap
    }
    
    # Return top-left and bottom-right corners of the overlapping area
    return(list(top_left=c(overlap_x1, overlap_y1),
                bottom_right=c(overlap_x2, overlap_y2)))
}


#################################################

# 1. READ RASTER DATA FROM 7 TXT FILES

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

# Crop sea areas (-999) to 0
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
rm(tenerife1, tenerife2, tenerife3, tenerife4, tenerife5, tenerife6, tenerife7)

# Put all 7 images together manually in Photoshop -> "tenerifecomposite.tif"
RESOLUTION=25
MAXIMO=3710.062
DEM=readTIFF("tenerifecomposite.tif")*MAXIMO
hist(DEM[DEM>0], breaks=800)


#################################################

# 2. RESCALE DEM TO FULL HD (1920 X 1080) AND CREATE SOLID MAP

f=2
fscale=0.347758887171561*f  # downsampling scaling (chosen to fit in Full HD)
RESOLUTION=RESOLUTION/fscale  # downsampling also reduces RESOLUTION
DIMY=round(nrow(DEM)*fscale)
DIMX=round(ncol(DEM)*fscale)
DEMrs=arrayresample(DEM, DIMX, DIMY)
# DEMrs[DEMrs<0]=0  # clip to 0 Lanczos artifacts

DIMYfullHD=1080*f
DIMXfullHD=1920*f
DEM=matrix(0, nrow=DIMYfullHD, ncol=DIMXfullHD)
DEM[(DIMYfullHD/2-DIMY/2):(DIMYfullHD/2+DIMY/2-1),
    (DIMXfullHD/2-DIMX/2):(DIMXfullHD/2+DIMX/2-1)]=DEMrs
rm(DEMrs)

writeTIFF(DEM/max(DEM), "tenerifecomposite_fullHD.tif",
          bits.per.sample=16, compression="LZW")

# Solid map
solid=DEM
solid[solid>0]=1
writeTIFF(solid, "solid.tif",
          bits.per.sample=16, compression="LZW")


#################################################

# 3. GENERATE STANDARD HILLSHADE

hillshade=hillshademap(DEM, dx=RESOLUTION, dlight=c(0, 1, 1))

# Save hillshade
writeTIFF(hillshade, "hillshade.tif", bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hillshade[nrow(hillshade):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=1)),
      asp=nrow(hillshade)/ncol(hillshade), axes=FALSE)


#################################################

# 4. DEM DISCRETIZATION

# Obtain blurred version of DEM
DEMblur=arrayblur(DEM)
DEMblur[DEM==0]=0  # preserve original coastline
DEMblur=DEMblur-min(DEMblur)
DEMblur=DEMblur/max(DEMblur)
writeTIFF(DEMblur, "tenerifecomposite_fullHD_blur.tif",
          bits.per.sample=16, compression="LZW")

# Slicing
NSLICES=20
DEMslice=DEMblur
DEMslice[DEMslice==0]=-1/NSLICES  # water areas
DEMslice=floor(DEMslice*NSLICES)/NSLICES
DEMslice[DEMslice==max(DEMslice)]=(NSLICES-1)/NSLICES

hist(DEMslice, breaks=200)
DEMslice=DEMslice-min(DEMslice)
DEMslice=DEMslice/max(DEMslice)
writeTIFF(DEMslice, paste0("tenerife_", NSLICES, "slices.tif"),
          bits.per.sample=16, compression="LZW")


#################################################

# 5. CALCULATE CONTOURS AND TANAKA HILLSHADE

# Contours
DEMcontour=contour(DEMslice, stroke=3)
writeTIFF(DEMcontour, paste0("tenerife_contours.tif"),
          bits.per.sample=16, compression="LZW")

# Tanaka contours: slight blur + new hillshade
DEMtanaka=arrayblur(DEMslice, radius=1)
writeTIFF(DEMtanaka, paste0("tenerife_tanaka.tif"),
          bits.per.sample=16, compression="LZW")

hillshadetanaka=hillshademap(DEMtanaka*MAXIMO, dx=RESOLUTION, dlight=c(0, 1, 1))
writeTIFF(hillshadetanaka, "hillshadetanaka.tif",
          bits.per.sample=16, compression="LZW")


#################################################

# 6. CALCULATE SHADOWS

shadows=shadowmap(DEM, dx=RESOLUTION, dlight=c(0, 30, 5))
writeTIFF(shadows, "shadowshard.tif", bits.per.sample=16, compression="LZW")

shadows=shadowmap(DEM, dx=RESOLUTION, dlight=c(0, 30, 5))
Gamma=6
writeTIFF(shadows^Gamma, "shadowsmixed.tif", bits.per.sample=16, compression="LZW")


#################################################

# 7. ANIM

# Read hillshade and lens shape
img1=readTIFF("background.tif")  # background hillshade
DIMY=nrow(img1)
DIMX=ncol(img1)

figure=readPNG("circlelarge.png")  # magnifying lens shape
dimy=nrow(figure)
dimx=ncol(figure)
RADIUSY=dimy/2
RADIUSX=dimx/2

# Read maps
mapas=c('DEM', 'contours', 'colours', 'shadows')
NMAPS=length(mapas)
img2=array(0, c(dim(img1), NMAPS))  # each of 4 partial maps
for (i in 1:NMAPS) img2[,,,i]=readTIFF(paste0(mapas[i], ".tif"))

x0=DIMX/2
y0=DIMY/2
R=320
NFRAMES=360  # total number of frames to calculate
for (frame in 0:(NFRAMES-1)) {
    THETA=frame*(2*pi)/NFRAMES
    imgout=img1  # reset to background
    for (i in 1:NMAPS) {
        theta=THETA+(i-1)*pi/2
        POSX=round(x0+R*cos(theta)-RADIUSX)  # last minute rounding
        POSY=round(y0+R*sin(theta)-RADIUSY)
        
        mask=img1*0
        SOLAPE=overlap(DIMX, DIMY, POSX, POSY, dimx, dimy)
        if (!is.null(SOLAPE)) {
            # Call again the same function swapping parameters
            solape=overlap(dimx, dimy, -POSX+2, -POSY+2, DIMX, DIMY)
            mask[SOLAPE[[1]][2]:SOLAPE[[2]][2],
                 SOLAPE[[1]][1]:SOLAPE[[2]][1],]=replicate(3,
                                       figure[solape[[1]][2]:solape[[2]][2],
                                              solape[[1]][1]:solape[[2]][1]])   
        }
        imgout = (1-mask)*imgout + mask*img2[,,,i]  # image masking  
    }
    name=paste0("imgout", ifelse(frame<10,"00",ifelse(frame<100,"0","")),
                frame, ".png")
    print(name)
    writePNG(imgout, name)  # save frame
}


# MP4 Video (MPEG-4 AVC/H.264):
# ffmpeg -loop 1 -framerate 24 -i imgout%03d.png -i interfrontstrangeclip.wav
# -t 96.872 -c:v libx264 -crf 20 -pix_fmt yuv420p xraymap.mp4

