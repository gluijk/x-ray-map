# Level curves representation
# www.overfitting.net
# https://www.overfitting.net/2017/12/visualizaciones-ad-hoc-de-sonido-con-r.html

library(terra)
library(tiff)
library(png)

hillshade=function(map, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # dx: map resolution (m)
    # Lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    
    DIMY=nrow(map)    
    DIMX=ncol(map)
    
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(map[1:(DIMY-2), 2:(DIMX-1)] - map[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(map[2:(DIMY-1), 1:(DIMX-2)] - map[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 'lost' borders
    hillshade=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshade[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshade[c(1,DIMY),]=hillshade[c(2,DIMY-1),]
    hillshade[,c(1,DIMX)]=hillshade[,c(2,DIMX-1)]
    
    return(hillshade^(1/gamma))
}


img=readTIFF("tenerife_expanded500px_blur20.tif")
# img=readTIFF("tenerife_expanded500px_blur20_1920.tif")
# img=readTIFF("guadarrama_blur5.tif")

img=img-min(img)
img=img/max(img)

DIMY=nrow(img)
DIMX=ncol(img)

# Slicing
NSLICES=20
img[img==0]=-1/NSLICES  # water areas
img=floor(img*NSLICES)/NSLICES
img[img==max(img)]=(NSLICES-1)/NSLICES

hist(img, breaks=200)
img=img-min(img)
img=img/max(img)
writeTIFF(img, paste0("tenerife_sliced",NSLICES,".tif"),
         bits.per.sample=16, compression="LZW")

# Pseudo 3D view
ALTTOTAL=420 # ALTTOTAL=DELTAZ*NSLICES
imgout=img*0  # array(0, c(DIMY*2, DIMX))

values=sort(unique(as.vector(img)))
NVALUES=length(values)  # NVALUES = NSLICES + 1
DELTAZ=round(ALTTOTAL/NSLICES)+1
type='ring' # solid'  # ring'
for (i in 2:NVALUES) {
    print(paste0(i, " de ", NVALUES))
    if (type=='ring') indices=which(img == values[i], arr.ind=TRUE)
        else indices=which(img >= values[i], arr.ind=TRUE)
    for (j in 1:DELTAZ) {
        indicesdest=indices
        indicesdest[,1]=indicesdest[,1]-DELTAZ*(i-1)+j
        imgout[indicesdest]=0.01
    }
    indicesdest=indices
    indicesdest[,1]=indicesdest[,1]-DELTAZ*(i-1)
    imgout[indicesdest]=values[i]
}

writePNG(imgout, paste0("stepsmap_", NSLICES, ".png"))


# Calculate outline map from solid map
outline=imgout*0
# 1 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=
    abs(imgout[1:(DIMY-2), 2:(DIMX-1)] -
            imgout[2:(DIMY-1), 2:(DIMX-1)]) +
    abs(imgout[2:(DIMY-1), 1:(DIMX-2)] -
            imgout[2:(DIMY-1), 2:(DIMX-1)])
# increase to 2 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
    outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
# increase to 3 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
    outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
    outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]

outline[outline!=0]=1
writePNG(outline, "mapoutline.png")

borders=which(outline==1, arr.ind=TRUE)
BORDERSLEN=nrow(borders)
aleat=round(runif(BORDERSLEN))
for (i in 1:BORDERSLEN) {
    if (aleat[i]) {
        y=borders[i,1]
        x=borders[i,2]
        outline[(y-1):(y+1), x]=1
        outline[y, (x-1):(x+1)]=1        
    }
}
borders=which(outline==1)  # update borders including noise
writePNG(outline, "mapoutline_noise.png")


imgout[borders]=0.01  # value for borders
imgout[imgout>0.01]=1
imgout[imgout==0]=0.5
imgout[imgout==0.01]=0
writePNG(imgout, paste0("stepsmap_", NSLICES, "_borders.png"))





# Hillshade from slicing
img=readTIFF("guadarrama_sliced10_blur2.tif")
hill=hillshade(img*1000, dx=25, dlight=c(0, 2, 3), gamma=1)
writeTIFF(hill, "hillshade.tif", bits.per.sample=16, compression="LZW")



###################################
# MINECRAFT STYLE TENERIFE


# Distorted points (source)
imgd=readTIFF("distorted.tif")
xu=c(1400, 1600, 1600, 1400)  # top-left, bottom-left, bottom-right, top-right
yu=c(1400, 1400, 1600, 1600)

# Undistorted points (destination)
# imgu=readTIFF("undistorted.tif")  # not used
xd=c(1500, 1600, 1500, 1400)  # top-left, bottom-left, bottom-right, top-right
yd=c(1400, 1500, 1600, 1500)

# NOTE: we swap the distorted and undistorted trapezoids because
# we want to model the transformation
# FROM CORRECTED coords (DST) -> TO UNCORRECTED coords (ORG)


# Solve 8 equations linear system: A * k = b -> k = inv(A) * b
A=matrix(nrow=8, ncol=8)
A[1,]=c(xd[1], yd[1], 1, 0,     0,     0, -xd[1]*xu[1], -yd[1]*xu[1])
A[2,]=c(0,     0,     0, xd[1], yd[1], 1, -xd[1]*yu[1], -yd[1]*yu[1])
A[3,]=c(xd[2], yd[2], 1, 0,     0,     0, -xd[2]*xu[2], -yd[2]*xu[2])
A[4,]=c(0,     0,     0, xd[2], yd[2], 1, -xd[2]*yu[2], -yd[2]*yu[2])
A[5,]=c(xd[3], yd[3], 1, 0,     0,     0, -xd[3]*xu[3], -yd[3]*xu[3])
A[6,]=c(0,     0,     0, xd[3], yd[3], 1, -xd[3]*yu[3], -yd[3]*yu[3])
A[7,]=c(xd[4], yd[4], 1, 0,     0,     0, -xd[4]*xu[4], -yd[4]*xu[4])
A[8,]=c(0,     0,     0, xd[4], yd[4], 1, -xd[4]*yu[4], -yd[4]*yu[4])

b=as.matrix(c(xu[1], yu[1], xu[2], yu[2], xu[3], yu[3], xu[4], yu[4]))

k=solve(A, b)  # equivalent to inv(A) * b = solve(A) %*% b

# Undo distortion function
undo.keystone = function(xd, yd, k) {
    xu=(k[1]*xd+k[2]*yd+k[3]) / (k[7]*xd+k[8]*yd+1)
    yu=(k[4]*xd+k[5]*yd+k[6]) / (k[7]*xd+k[8]*yd+1)
    return(c(xu, yu))  # return pair (xu, yu)
}

# Check
for (i in 1:4) {
    print(undo.keystone(xd[i], yd[i], k))
}


# Plot trapezoids
plot(c(xd, xd[1]), c(yd, yd[1]), type='l', col='red', asp=1,
     xlab='X', ylab='Y', xlim=c(0, 4000), ylim=c(6000, 0))
lines(c(xu, xu[1]), c(yu, yu[1]), type='l', col='blue')
for (i in 1:4) {
    lines(c(xd[i], xu[i]), c(yd[i], yu[i]), type='l', lty=3, col='darkgray')
}
abline(h=c(0,6000), v=c(0,4000))


# Correct keystone distortion
DIMXd=ncol(imgd)
DIMYd=nrow(imgd)

EDGE=1000  # additional edge (Y axis)
imgc=array(0, dim=c(DIMYd+EDGE, DIMXd))  # imgc=imgd*0

DIMXc=ncol(imgc)
DIMYc=nrow(imgc)
for (x in 1:DIMXc) {
    for (y in 1:DIMYc) {
        xuyu=round(undo.keystone(x, y+EDGE*0.2, k))  # add bottom
        if (xuyu[1]>=1 & xuyu[1]<=DIMXd & xuyu[2]>=1 & xuyu[2]<=DIMYd)
            imgc[y, x]=imgd[xuyu[2], xuyu[1]]  # nearest neighbour interp
    }
}
writeTIFF(imgc, "corrected.tif", bits.per.sample=16)


img=readTIFF("corrected.tif")

img=img-min(img)
img=img/max(img)

DIMY=nrow(img)
DIMX=ncol(img)

# Slicing
NSLICES=20
img[img==0]=-1/NSLICES  # water areas
img=floor(img*NSLICES)/NSLICES
img[img==max(img)]=(NSLICES-1)/NSLICES

hist(img, breaks=200)
img=img-min(img)
img=img/max(img)
writeTIFF(img, paste0("tenerife_sliced",NSLICES,".tif"),
          bits.per.sample=16, compression="LZW")

# Pseudo 3D view
ALTTOTAL=150 # ALTTOTAL=DELTAZ*NSLICES
imgout=img*0  # array(0, c(DIMY*2, DIMX))

values=sort(unique(as.vector(img)))
NVALUES=length(values)  # NVALUES = NSLICES + 1
DELTAZ=round(ALTTOTAL/NSLICES)+1
type='ring' # solid'  # ring'
for (i in 2:NVALUES) {
    print(paste0(i, " de ", NVALUES))
    if (type=='ring') indices=which(img == values[i], arr.ind=TRUE)
    else indices=which(img >= values[i], arr.ind=TRUE)
    for (j in 1:DELTAZ) {
        indicesdest=indices
        indicesdest[,1]=indicesdest[,1]-DELTAZ*(i-1)+j
        imgout[indicesdest]=0.01
    }
    indicesdest=indices
    indicesdest[,1]=indicesdest[,1]-DELTAZ*(i-1)
    imgout[indicesdest]=values[i]
}

writePNG(imgout, paste0("stepsmap_", NSLICES, ".png"))



# Calculate outline map from solid map
outline=imgout*0
# 1 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=
    abs(imgout[1:(DIMY-2), 2:(DIMX-1)] -
            imgout[2:(DIMY-1), 2:(DIMX-1)]) +
    abs(imgout[2:(DIMY-1), 1:(DIMX-2)] -
            imgout[2:(DIMY-1), 2:(DIMX-1)])
# increase to 2 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
    outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
# increase to 3 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
    outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
    outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]

outline[outline!=0]=1
writePNG(outline, "mapoutline.png")

borders=which(outline==1, arr.ind=TRUE)
BORDERSLEN=nrow(borders)
aleat=round(runif(BORDERSLEN))
for (i in 1:BORDERSLEN) {
    if (aleat[i]) {
        y=borders[i,1]
        x=borders[i,2]
        outline[(y-1):(y+1), x]=1
        outline[y, (x-1):(x+1)]=1        
    }
}
borders=which(outline==1)  # update borders including noise
writePNG(outline, "mapoutline_noise.png")


imgout[borders]=0.01  # value for borders
imgout[imgout>0.01]=1
imgout[imgout==0]=0.5
imgout[imgout==0.01]=0
writePNG(imgout, paste0("stepsmap_", NSLICES, "_borders.png"))


