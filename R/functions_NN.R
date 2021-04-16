SampleNorm <- function(train_count,
                       test_count = NULL,
                       nfeature = 5000) {

  ### sample normalization
  cmax <- apply(train_count, 2, max)
  cmin <- apply(train_count, 2, min)
  tmp1 <- sweep(train_count, 2, cmin, "-")
  tmp2 <- sweep(tmp1, 2, cmax-cmin, "/")

  ### let's try some feature selection here
  oosd <- apply(tmp2, 1, sd)
  tmp3 <- tmp2[order(oosd, decreasing = TRUE)[1:nfeature], ]
  train_count_out <- t(tmp3)
  rm(tmp1, tmp2, tmp3)

  if(!is.null(test_count)) {
    ### sample normalization
    cmax <- apply(test_count, 2, max)
    cmin <- apply(test_count, 2, min)
    tmp1 <- sweep(test_count, 2, cmin, "-")
    tmp2 <- sweep(tmp1, 2, cmax-cmin, "/")

    ### let's try some feature selection here
    cname <- intersect(colnames(train_count_out), rownames(tmp2))
    tmp3 <- tmp2[cname, ]
    test_count_out <- t(tmp3)
    rm(tmp1, tmp2, tmp3)
    train_count_out = train_count_out[,cname]
    return(list(train_count_out = train_count_out,
                test_count_out = test_count_out))
  } else {
    return(list(train_count_out = train_count_out))
  }
}

BigNN <- function(train_count, train_label,
                  lossname = 'categorical_crossentropy',
                  last_act = "softmax",
                  nClass = 8,
                  verbose = TRUE) {
  require(keras)

  ### Build a big neural network:
  model_big <- keras_model_sequential() %>%
    layer_dense(units = 256, activation = "relu", input_shape = ncol(train_count)) %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = nClass, activation = last_act) %>%
    compile(
      loss = lossname,
      optimizer = optimizer_rmsprop(),
      metrics = c('accuracy')
    )
  model_big %>%
    fit(
      x = train_count,
      y = train_label,
      epochs = 15,
      batch_size = 256,
      validation_split =0.2,
      # callbacks = list(
      #   callback_early_stopping(patience = 10),
      #   callback_reduce_lr_on_plateau(factor = 0.05)
      # ),
      verbose = verbose,
      return_best_model=TRUE
    )

  return(model_big)
}


MediumNN <- function(train_count, train_label,
                     lossname = 'categorical_crossentropy',
                     last_act = "softmax",
                     nClass = 8,
                     verbose = TRUE) {
  require(keras)

  ### Build a big neural network:
  model_medium <- keras_model_sequential() %>%
    layer_dense(units = 128, activation = "relu", input_shape = ncol(train_count)) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = nClass, activation = last_act) %>%
    compile(
      loss = lossname,
      optimizer = optimizer_rmsprop(),
      metrics = c('accuracy')
    )
  model_medium %>%
    fit(
      x = train_count,
      y = train_label,
      epochs = 15,
      batch_size = 256,
      validation_split =0.2,
      callbacks = list(
        callback_early_stopping(patience = 10),
        callback_reduce_lr_on_plateau(factor = 0.05)
      ),
      verbose = verbose,
      return_best_model=TRUE
    )

  return(model_medium)
}


SmallNN <- function(train_count, train_label,
                    lossname = 'categorical_crossentropy',
                    last_act = "softmax",
                    nClass = 8,
                    verbose = TRUE) {
  require(keras)

  ### Build a big neural network:
  model_small <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu", input_shape = ncol(train_count)) %>%
    layer_dense(units = nClass, activation = last_act) %>%
    compile(
      loss = lossname,
      optimizer = optimizer_rmsprop(),
      metrics = c('accuracy')
    )
  model_small %>%
    fit(
      x = train_count,
      y = train_label,
      epochs = 15,
      batch_size = 256,
      validation_split =0.2,
      verbose = verbose,
      return_best_model=TRUE
    )

  return(model_small)
}

