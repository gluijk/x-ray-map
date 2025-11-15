# Beautiful map of El Hierro (the Canary Islands)
# www.overfitting.net
# https://www.overfitting.net/2024/10/radiografia-de-tenerife-con-r.html

library(terra)  # build blur and resample functions
library(tiff)  # save 16-bit TIFF's


# Blur
# https://stackoverflow.com/questions/70429190/how-can-i-perform-neighborhood-analysis-in-terra-or-raster-and-keep-the-same-na
arrayblur=function(img, radius=11) {
    # radius: radius of the circular averaging window
    require(png)
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


# Hillshade calculation
hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction. It can be defined in three ways:
    #   a) 1D value indicating light source Azimuth in degrees (0-360)
    #      (0=North, 90=East, 180=South, 270=West)
    #   b) 2D vector indicating light source (X,Y) coordinates
    #   c) 3D vector indicating light source (X,Y,Z) coordinates:
    #      (X=South, Y=East, Z=Up)
    #      dlight=c(0, 2, 3)  # sunrise
    #      dlight=c(0, 0, 1)  # midday
    #      dlight=c(0,-2, 3)  # sunset
    #   NOTE: both in a) and b) a 45º Elevation angle is applied
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    # Deal with lighting direction
    if (length(dlight)==1) dlight=c(-cos(dlight*pi/180), sin(dlight*pi/180))
    if (length(dlight)==2) dlight=c(dlight, (dlight[1]^2+dlight[2]^2)^0.5)
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


#################################################

# 1. PROCESS GEOTIFF DATA AND CONVERT FROM DEGREES TO M

# Read 2m resolution GeoTIFF files
hierro=mosaic(rast("MDT02-WGS84-1105-2-COB2.tif"),
              rast("MDT02-WGS84-1105-3-COB2.tif"),
              rast("MDT02-WGS84-1105-4-COB2.tif"),
              rast("MDT02-WGS84-1108-2-1-COB2.tif"),
              fun='mean')
hierro
plot(hierro)


# REPROJECT to m (ej. REGCAN95 / UTM 28N for Canary Islands)
hierro=project(hierro, "EPSG:4083")  # takes long...
hierro
plot(hierro)


#################################################

# 2. CROP DEM, RESAMPLE TO FULL HD (1920 X 1080) AND CREATE SOLID MAP

# EXTEND a bit (+2km) and CROP raster to area of interest
e <- ext(hierro)
e_new <- ext(e$xmin, e$xmax + 2000, e$ymin, e$ymax)
hierro <- extend(hierro, e_new)

cropdef=ext(187000, 219000, 3059000, 3086000)
hierrocrop=crop(x=hierro, y=cropdef)
hierrocrop
plot(hierrocrop)


# RESAMPLE to Full HD
DIMX=1920*2
DIMY=round(DIMX*nrow(hierrocrop)/ncol(hierrocrop))

hierrocroprs=rast(nrows=DIMY, ncols=DIMX, extent=ext(hierrocrop))
hierrocroprs=resample(x=hierrocrop, y=hierrocroprs, method='bilinear', threads=TRUE)
hierrocroprs
plot(hierrocroprs)

RESOLUTION=res(hierrocroprs)[1]
# Map scale
print(paste0("5km over the map of width=", DIMX, " pixels correspond to ",
             round(5000/RESOLUTION), " pixels"))
MAXIMO=as.integer(global(hierrocroprs, "max", na.rm=TRUE))


# Convert to matrix and save as TIFF
DEM=as.matrix(hierrocroprs, wide=TRUE)
DEM[is.na(DEM)]=0
DEM[DEM<3]=0  # delete some residuals along the coast

v=DEM[DEM>0]
hist(v,
     breaks=100, xlim=c(0, MAXIMO),
     main=paste0("Distr. altitudes El Hierro"),
     xlab=paste0("min / mediana / media / max = ",
                 round(min(v)*0), "m / ",
                 round(median(v)), "m / ",
                 round(mean(v)), "m / ",
                 round(max(v)), "m"))
abline(v=median(v), col='red', lty='dashed', lwd=2)
abline(v=mean(v), col='red', lty='dashed', lwd=2)

writeTIFF(DEM/max(DEM), "dem.tif", bits.per.sample=16, compression='LZW')


# Solid map
solid=DEM
solid[solid>0]=1
writeTIFF(solid, "solid.tif", bits.per.sample=16, compression="LZW")


#################################################

# 3. GENERATE STANDARD HILLSHADE

# dlight(X=South, Y=East, Z=Up)
hillshade=hillshademap(DEM, dx=RESOLUTION*1.8, dlight=c(0, -2, 2))+
    hillshademap(DEM, dx=RESOLUTION, dlight=c(0, --4, 2))*0.2*0
hillshade=hillshade/max(hillshade)

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
writeTIFF(DEMblur, "hierro_blur.tif",
          bits.per.sample=16, compression="LZW")

# Slicing
NSLICES=12
DEMslice=DEMblur
DEMslice[DEMslice==0]=-1/NSLICES  # water areas
DEMslice=floor(DEMslice*NSLICES)/NSLICES
DEMslice[DEMslice==max(DEMslice)]=(NSLICES-1)/NSLICES

hist(DEMslice, breaks=200)
DEMslice=DEMslice-min(DEMslice)
DEMslice=DEMslice/max(DEMslice)
writeTIFF(DEMslice, paste0("hierro_", NSLICES, "slices.tif"),
          bits.per.sample=16, compression="LZW")


#################################################

# 5. CALCULATE CONTOURS AND TANAKA HILLSHADE

# Contours
DEMcontour=contour(DEMslice, stroke=3)
writeTIFF(DEMcontour, paste0("hierro_contours.tif"),
          bits.per.sample=16, compression="LZW")

# Tanaka contours: slight blur + new hillshade
DEMtanaka=arrayblur(DEMslice, radius=2)
writeTIFF(DEMtanaka, paste0("hierro_tanaka.tif"),
          bits.per.sample=16, compression="LZW")

hillshadetanaka=hillshademap(DEMtanaka*MAXIMO, dx=RESOLUTION, dlight=c(0, 1, 1))
writeTIFF(hillshadetanaka, "hillshadetanaka.tif",
          bits.per.sample=16, compression="LZW")


#################################################

# 6. CALCULATE SHADOWS

shadows=shadowmap(DEM, dx=RESOLUTION, dlight=c(0, -30, 5))
writeTIFF(shadows, "shadowshard.tif", bits.per.sample=16, compression="LZW")

shadows=shadowmap(DEM, dx=RESOLUTION, dlight=c(0, 15, 5))
Gamma=6
writeTIFF(shadows^Gamma, "shadowsmixed.tif", bits.per.sample=16, compression="LZW")
