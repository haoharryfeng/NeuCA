GetSensitivity <- function(mycount, onelabel) {
  sensvec <- rep(NA, nrow(mycount))
  TP <- rowSums(mycount[, onelabel == 1]>0)
  AllP <- rowSums(mycount>0)
  sensvec <- TP/AllP
  return(sensvec)
}

GetAllMetric <- function(mycount, truelabels) {
  message("Calculating gene-specific sensitivities...")

  uniqCT <- unique(truelabels)
  allmet <- matrix(NA, nrow(mycount), length(uniqCT))

  for(i in seq_along(uniqCT)) {
    onelabel <- rep(0, length(truelabels))
    onelabel[truelabels == uniqCT[i]] <- 1

    allmet[, i] <- GetSensitivity(mycount, onelabel)
  }
  onelabel <- rep(0, length(truelabels))

  colnames(allmet) <- paste0(uniqCT, "_Sensitivity")
  rownames(allmet) <- rownames(mycount)
  return(allmet)
}

SelectMarker <- function(mytrain,
                         trainlabel,
                         gm_threshold = 20,
                         resout) {
  message("Selecting best markers...")

  uniqCT <- unique(trainlabel)

  ### Get markers based on sensitivity
  rownames(resout) <- rownames(mytrain)
  gmidx <- which(rowSums(mytrain) > gm_threshold)
  allsens <- resout[gmidx, ]

  bestmarker <- list()
  allsensvec <- rep(0, length(uniqCT))
  for(i in seq_along(uniqCT)) {
    tmp0 <- allsens[order(allsens[,i], decreasing = TRUE),]
    tmp1 <- tmp0[tmp0[,i]>0.95,]
    if(is.matrix(tmp1) & length(tmp1)>0) {
      bestmarker[[i]] <- rownames(tmp1)
      allsensvec[i] <- mean(tmp1[,i])
    } else {
      bestmarker[[i]] <- rownames(tmp0)[1:5]
      allsensvec[i] <- mean(tmp0[1:5,1])
    }
  }
  names(bestmarker) <- uniqCT
  names(allsensvec) <- uniqCT
  allsensvec2 <- allsensvec[order(allsensvec)]
  bestmarker2 <- bestmarker[order(allsensvec)]

  return(list(allsensvec = allsensvec2,
              bestmarker = bestmarker2))
}

assignCT <- function(mycount,
                     allsensvec2,
                     bestmarker2,
                     csthreshold = 3) {
  message("Assignning cell types based on marker genes...")

  estlabels <- rep(NA, ncol(mycount))

  for(i in seq_along(allsensvec2)) {
    nidx <- match(bestmarker2[[i]], rownames(mycount))
    onecount <- mycount[na.omit(nidx), ]
    if(!is.null(nrow(onecount))) {
      cs <- colSums(onecount)
    } else {
      cs <- onecount
    }
    estlabels[which(cs>=csthreshold)] <- names(allsensvec2)[i]
  }
  return(estlabels)
}

GetFeature <- function(LeftRef, RightRef, nMark) {
  require(limma)

  numlab <- c(rep(1, ncol(LeftRef)), rep(0, ncol(RightRef)))
  allRef <- cbind(LeftRef, RightRef)
  allRef2 <- allRef[which(rowMeans(allRef)>0),]
  fit <- lmFit(allRef2, numlab)
  fit <- eBayes(fit)
  fit2 <- topTable(fit, number = nrow(allRef2))

  return(rownames(fit2)[1:nMark])
}

GetTree <- function(uniqlab, treestr) {
  treelabels <- matrix(NA, nrow(treestr), ncol(treestr))
  treelabels[treestr<0] <- uniqlab[-treestr[treestr<0]]
  treelabels[treestr>0] <- paste0("Step", treestr[treestr>0])

  tmplab <- treelabels
  for(ii in 1:nrow(treestr)) {
    for(jj in 1:2) {
      if(treestr[ii,jj]>0) {
        tmplab[ii,jj] <- paste0(tmplab[treestr[ii,jj],], collapse = " ")
      }
    }
  }

  treeCellTypes <- list()
  for(ii in 1:nrow(treestr)) {
    treeCellTypes[[ii]] <- list()
    for(jj in 1:2) {
      tmptt <- unlist(strsplit(tmplab[ii,jj], split = " "))
      names(tmptt) <- NULL
      treeCellTypes[[ii]][[jj]] <- tmptt
    }
  }

  return(list(treelabels = treelabels[rev(1:nrow(treelabels)), ],
              treeCellTypes = rev(treeCellTypes)))
}

trainOneNN <- function(LeftRef, RightRef,
                       modeltype = "medium",
                       verbose = FALSE) {
  require(e1071)
  numlab <- c(rep(1, ncol(LeftRef)), rep(0, ncol(RightRef)))
  inputcount <- as.matrix(t(cbind(LeftRef, RightRef)))

  train_y <- to_categorical(numlab,2)

  if (modeltype == "big") {
    outmodel <- BigNN(inputcount, train_y, lossname = 'binary_crossentropy', last_act = "sigmoid", nClass = 2, verbose = verbose)
  } else if (modeltype == "medium") {
    outmodel <- MediumNN(inputcount, train_y, lossname = 'binary_crossentropy', last_act = "sigmoid", nClass = 2, verbose = verbose)
  } else if (modeltype == "small") {
    outmodel <- SmallNN(inputcount, train_y, lossname = 'binary_crossentropy', last_act = "sigmoid", nClass = 2, verbose = verbose)
  }

  return(outmodel)
}

train_SCH <- function(train_count, train_lab,
                      nMark = 1000,
                      modeltype = modeltype,
                      verbose = TRUE) {

  allfeatures <- list()

  ### get cell type mean profile
  uniqlab <- unique(train_lab)
  grpmean <- matrix(NA, nrow(train_count), length(uniqlab))
  for(i in seq_along(uniqlab)) {
    grpmean[,i] <- rowMeans(train_count[,which(train_lab == uniqlab[i])])
  }
  colnames(grpmean) <- uniqlab
  # grpmean2 <- grpmean[,]
  grpmean <- na.omit(grpmean)

  ### hierarchical clustering
  dd <- dist(t(grpmean))
  hh <- hclust(dd)
  Htree <- GetTree(uniqlab = hh$labels,
                   treestr = hh$merge)
  treelabels <- Htree$treelabels
  treeCellTypes <- Htree$treeCellTypes

  ### get hierarchical NN tree
  HNNtrees <- list()
  for(iter in 1:nrow(treelabels)) {
    message(paste0("Training hierarchical NN: Step", iter))
    LeftRef <- train_count[, train_lab %in% treeCellTypes[[iter]][[1]]]
    RightRef <- train_count[, train_lab %in% treeCellTypes[[iter]][[2]]]
    thisfeature <- GetFeature(LeftRef, RightRef, nMark = nMark)
    allfeatures[[iter]] <- thisfeature
    matidx1 <- match(thisfeature, rownames(LeftRef))
    matidx2 <- match(thisfeature, rownames(RightRef))
    LeftRef2 <- LeftRef[matidx1, ]
    RightRef2 <- RightRef[matidx2, ]
    HNNtrees[[iter]] <- trainOneNN(LeftRef2, RightRef2,
                                   modeltype = modeltype,
                                   verbose = verbose)
  }

  return(list(treelabels = treelabels,
              treeCellTypes = treeCellTypes,
              HNNtrees = HNNtrees,
              allfeatures = allfeatures))
}


HNN_assign <- function(input_count, estlabels, train_SCH_out) {

  treelabels = train_SCH_out$treelabels
  treeCellTypes = train_SCH_out$treeCellTypes
  HNNtrees = train_SCH_out$HNNtrees
  allfeatures = train_SCH_out$allfeatures
  rownames(treelabels) <- paste0("Step", nrow(treelabels):1)

  message("Assign labels through the hierarchical Neural Network tree:")
  message(paste0("A total of ", nrow(treelabels), "steps!"))
  pb <- txtProgressBar(min = 0, max = nrow(treelabels))
  for(iter in 1:nrow(treelabels)) {
    setTxtProgressBar(pb, iter)
    input_count2 <- input_count[allfeatures[[iter]], ]
    if(iter == 1) {
      int_indx <- which(is.na(estlabels))
    } else {
      int_indx <- which(estlabels == rownames(treelabels)[iter])
    }
    tmpcell <- as.matrix(input_count2[, int_indx])
    if(!is.null(ncol(tmpcell)) & is.matrix(tmpcell) & length(int_indx)>0) {
      tmplab <- rep(NA, ncol(tmpcell))
      thismodel <- HNNtrees[[iter]]
      pred_y <- predict(thismodel, x = as.matrix(t(tmpcell)),
                        batch_size = 64,
                        verbose = 0)

      pred <- predict(thismodel, t(as.matrix(tmpcell)), decision.value = TRUE)
      newpred <- apply(pred, 1, which.max) - 1
      tmplab[which(newpred == 1)] <- treelabels[iter, 1]
      tmplab[which(newpred == 0)] <- treelabels[iter, 2]
      estlabels[int_indx] <- tmplab
    } else if(is.vector(tmpcell)) {
      tmplab <- NA
      thismodel <- HNNtrees[[iter]]
      tmpcell2 <- matrix(tmpcell,ncol = 1)
      rownames(tmpcell2) <- allfeatures[[iter]]
      pred <- predict(thismodel, t(tmpcell2), decision.value = TRUE)
      tmplab[which(pred == 1)] <- treelabels[iter, 1]
      tmplab[which(pred == 0)] <- treelabels[iter, 2]
      estlabels[int_indx] <- tmplab
    }
  }
  close(pb)

  return(estlabels)
}

SCHwrapper <- function(train_count, train_label,
                       test_count,
                       nMark = 500,
                       modeltype = "medium",
                       CellTypeMark_thres = 3,
                       GeneMarkerExpr_thres = 10,
                       verbose = TRUE){
  resout <- GetAllMetric(mycount = train_count, train_label)
  MarkOut <- SelectMarker(mytrain = train_count,
                          trainlabel = train_label,
                          gm_threshold = GeneMarkerExpr_thres,
                          resout = resout)
  estlabels <- assignCT(test_count, MarkOut$allsensvec,
                        MarkOut$bestmarker, csthreshold = CellTypeMark_thres)

  train_count_normed <- train_count/max(train_count)
  test_count_normed <- test_count/max(test_count)
  oosd <- apply(train_count_normed, 1, sd)
  train_count_normed <- train_count_normed[order(oosd, decreasing = TRUE)[1:5000], ]
  matidx <- match(rownames(train_count_normed), rownames(test_count_normed))
  test_count_normed <- test_count_normed[matidx, ]

  train_SCH_out <- train_SCH(train_count = train_count_normed,
                             train_lab = train_label,
                             nMark = nMark,
                             modeltype = modeltype,
                             verbose = verbose)
  SCHout <- HNN_assign(test_count_normed,
                       estlabels,
                       train_SCH_out)
  return(SCHout)
}

