###Main wrapper function NeuCA and some additional utility functions
NeuCA = function(train, test, model.size = "big", verbose = FALSE){
    message("Working on scRNA-seq data cell label training and testing:")
    if (missing(train))
        stop("Training scRNA-seq data must be provided.")
    if (missing(test))
        stop("Testing scRNA-seq data must be provided.")
    if (!is(train, "SingleCellExperiment"))
        stop("Training scRNA-seq data must be a SingleCellExperiment
            object, see package vignette.")
    if (!is(test, "SingleCellExperiment"))
        stop("Testing scRNA-seq data must be a SingleCellExperiment
            object, see package vignette.")
    if (ncol(assay(train)) != nrow(colData(train)))
        stop("The number of cells (columns) in the SCE count matrix must match
             the number of rows in SCE colData labels.")
    if (!model.size %in% c("big", "medium", "small"))
        stop("Model size must be one of the followings: 'big', 'medium', 'small'")
    if (!is.logical(verbose))
        stop("verbose must be binary Boolean")


### Data normalization
normedData <- SampleNorm(train_count = assay(train),
                            test_count = assay(test),
                             nfeature = 5000)

### Determine if cor>0.95 exist
cd = cor.det(t(normedData$train_count_out), colData(train)[,1])

CTname = unique(colData(train)[,1])

tr.lab = matrix(0,
                nrow = nrow(normedData$train_count_out),
                ncol = length(CTname)
                )
for(i in seq_along(CTname)){
        idx = which(colData(train)[,1] == CTname[i])
        tr.lab[idx, i] = 1
    }




###1. route of NN (cd = 1)
    if (cd == 1){
        message("Based on correlation values, direct neural network IS adopted")
        if (model.size == "big") {
            message("Neural network: big")
            outmodel <- BigNN(train_count = normedData$train_count_out,
                              train_label = tr.lab,
                              lossname = 'binary_crossentropy',
                              last_act = "sigmoid",
                              nClass = length(CTname),
                              verbose = verbose)
        } else if (model.size == "medium") {
            message("Neural network: medium")
            outmodel <- MediumNN(train_count = normedData$train_count_out,
                                 train_label = tr.lab,
                                 lossname = 'binary_crossentropy',
                                 last_act = "sigmoid",
                                 nClass = length(CTname),
                                 verbose = verbose)
        } else if (model.size == "small") {
            message("Neural network: small")
            outmodel <- SmallNN(train_count = normedData$train_count_out,
                                train_label = tr.lab,
                                lossname = 'binary_crossentropy',
                                last_act = "sigmoid",
                                nClass = length(CTname),
                                verbose = verbose)
        }


        preres <- predict(outmodel, x = as.matrix(normedData$test_count_out),
                          batch_size = 256,
                          verbose = verbose)

      prenum = apply(preres, 1, findm)
      predict.label = CTname[prenum]

    } else if (cd == 2){
        ###2. route of mhNN (cd = 2)
        message("Based on correlation values,
                marker-guided hierarchical neural network IS adopted")
        message("Marker-guided hierarchical neural network: ", model.size)
        predict.label <- SCHwrapper(train_count = t(normedData$train_count_out),
                           train_label = colData(train)[,1],
                           test_count = t(as.matrix(normedData$test_count_out)),
                           nMark = 500,
                           modeltype = model.size,
                           CellTypeMark_thres = 3,
                           GeneMarkerExpr_thres = 20,
                           verbose = verbose)
        }
    return(predict.label)
}


#this function determine the cell-type mean expression profile correlation
#return 1 if all cor < 0.95, return 2 of at least 1 cor >= 0.95
cor.det = function(dat, lb){
    #dat is a SingleCellExperiment training object
    ct = dat
    lb = lb
    lname = unique(lb)
    mp.ge = matrix(0, nrow = nrow(ct), ncol = length(lname))
    for(i in seq_along(lname)){
        #extract all cells GE, for each cell type category
        cge = ct[, which(lb == lname[i])]
        #average the profiles
        tmp = rowMeans(cge, na.rm = TRUE)
        #store them in mp.ge
        mp.ge[,i] = tmp
    }
    ac = cor(mp.ge)
    diag(ac) = 0
    #return 1 if all cor < 0.95, return 2 of at least 1 cor >= 0.95
    cd = ifelse(sum(ac<0.95) == length(ac), 1, 2)
    return(cd)
}

#####this function finds the index of elements that max a vector
findm = function(a){
    tmp = which(a == max(a))
    return(tmp)
}
