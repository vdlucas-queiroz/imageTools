library(terra)

# Funções relacionadas com o CMAP

#' Create Transition Matrix
#'
#' This function creates a transition matrix from multiple input files and writes the result to an output file.
#'
#' @param file_paths A character vector containing the file paths of the input transition matrices.
#' @param output A character string specifying the file path for the output transition matrix.
#' @return A data frame representing the combined transition matrix.
#' @importFrom utils read.table write.table
#' @export
#'
#' @examples
#' # Assuming example transition matrices are saved in "trans1.txt", "trans2.txt"
#' write.table(matrix(c(0.1, 0.9, 0.3, 0.7), nrow = 2), "trans1.txt", sep = ";", col.names = FALSE, row.names = FALSE)
#' write.table(matrix(c(0.8, 0.2, 0.6, 0.4), nrow = 2), "trans2.txt", sep = ";", col.names = FALSE, row.names = FALSE)
#' file_paths <- c("trans1.txt", "trans2.txt")
#' output <- "combined_transitions.txt"
#' transition_matrix(file_paths, output)
transition_matrix <- function(file_paths, output){

  num_transitions <- length(file_paths)
  num_classes_per_transition <- numeric(num_transitions + 1)
  transition_classes <- vector("list", num_transitions + 1)

  for (i in 1:num_transitions) {
    transition_matrix <- read.table(file_paths[i], sep=";", header=FALSE)

    num_classes_per_transition[i] <- dim(transition_matrix)[1]
    transition_classes[[i]] <- 1:num_classes_per_transition[i]

    if (i == num_transitions) {
      num_classes_per_transition[i + 1] <- dim(transition_matrix)[2]
      transition_classes[[i + 1]] <- 1:num_classes_per_transition[i + 1]
    }
  }

  all_combinations <- expand.grid(transition_classes)
  transition_weights <- numeric(nrow(all_combinations))

  for (i in 1:num_transitions) {
    transition_matrix <- read.table(file_paths[i], sep=";", header=FALSE)

    for (row in 1:nrow(transition_matrix)) {
      for (col in 1:ncol(transition_matrix)) {
        index <- which(all_combinations[, i] == row & all_combinations[, i + 1] == col)
        transition_weights[index] <- transition_matrix[row, col]
      }
    }

    all_combinations <- cbind(all_combinations, transition_weights)
  }

  if(num_transitions > 1) {

    # Calcula o número de colunas relacionadas a pesos na tabela 'tab'
    num_weight_cols <- num_transitions

    # Seleciona apenas as colunas de peso em 'tab'
    weight_cols <- (ncol(all_combinations) - num_weight_cols + 1):ncol(all_combinations)

    # Calcula o produto dos pesos em cada linha
    final_weights <- apply(all_combinations[, weight_cols], 1, prod)

    all_combinations <- all_combinations[, -weight_cols]

    all_combinations <- cbind(all_combinations, final_weights)
  }

  # Escreve a tabela final em um arquivo
  write.table(all_combinations, output, sep=";")
  print('TXT de transições criado...')
  return(all_combinations)
}

#sintaxe----------
# output_C3 <-"G:/Meu Drive/INPE/projeto_dissertacao/7_CMAP/C3/pesos_C3.txt"
#
# file_paths_C3<-c("G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2005-2006.txt",
#                  "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2006-2007.txt",
#                  "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2007-2008.txt",
#                  "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2008-2009.txt",
#                  "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2009-2010.txt")

# transition_matrix_C3 = transition_matrix(file_paths_C3,output_C3)



#-------------------------------------------------------------------------------------
# Function: CMAP_classifier
# Description: This function performs time-series Compound Maximum a Posteriori
# classification of raster images based on class likelihoods and trajectory weights.
# It computes the classification for each time point and outputs the result as raster files.

# Parameters:
#   class_likelihood_paths: A character vector containing the file paths to raster images for each date. Each raster file represents the likelihood of classes for a specific date.
#   trajectory_matrix: A dataframe representing the trajectories and their corresponding weights for classification.
#   output_path: A character string indicating the directory where the output raster files will be saved.

# Usage:
#   CMAP_classifier(class_likelihood_paths, trajectory_matrix, output_path)

# Example:
#   class_likelihood_paths <- c("path/to/2005_probs.tif", "path/to/2006_probs.tif")
#   trajectory_matrix <- read.table("path/to/trajectory_weights.txt", sep = ";", header = TRUE)
#   output_path <- "path/to/output/directory/"
#   CMAP_classifier(class_likelihood_paths, trajectory_matrix, output_path)

# Notes:
# - The function assumes that each raster file contains the likelihood for different classes at a given date.
# - The trajectory matrix should contain weights for different classification trajectories.
# - algorithm classifies all the images in the time series at the same time.
# - The output raster files will represent the classification result for each date and will be saved in the specified output directory.
# - The order of classes in the trajectory matrix should match the class IDs in the raster files.

# Author: [Mariane Souza Reis]
# Adaptation and optimization: [Vinícius D'Lucas Bezerra e Queiroz]
# Date: [30/01/2024]

#-------------------------------CMAP CLASSIFIER-------------------------------------


#' Compound Maximum A Posteriori (CMAP) Classifier
#'
#' This function implements the CMAP classifier, which uses class likelihoods and trajectory weights to classify remote sensing images.
#'
#' @param class_likelihood_paths A character vector containing the file paths of the class likelihood raster images for each date.
#' @param traj_weights A data frame containing the trajectory weights and the sequences of class transitions.
#' @param output_path A character string specifying the directory where the output classified images will be saved.
#' @param masking An optional character vector containing file paths of masking images for each date. Default is NULL.
#' @return A matrix representing the classification results for each date.
#' @importFrom terra rast values nlyr setValues writeRaster
#' @importFrom raster raster brick
#' @export
#'
CMAP_classifier <- function(class_likelihood_paths,traj_weights,output_path, masking = NULL ){


  n_dates <- length(class_likelihood_paths)
  n_classes <- unname(sapply(class_likelihood_paths, function(x) nlyr(rast(x))))

  # Pre-allocating the 'data' df with the correct size

  total_cols <- sum(n_classes)
  num_rows <- nrow(rast(class_likelihood_paths[1]))*ncol(rast(class_likelihood_paths[1]))

  data <- as.data.frame(matrix(nrow = num_rows, ncol = total_cols))

  current_col <- 1

  # Turning each likelihood image into a column of a df

  for (i in 1:n_dates){
    temp <- as.data.frame(values(rast(class_likelihood_paths[i])))

    # Determine the range of columns for 'temp' in 'data'
    next_col <- current_col + ncol(temp) - 1

    # Assign 'temp' to the range of columns in 'data'
    data[, current_col:next_col] <- temp

    # Update the index of the current column for the next iteration
    current_col <- next_col + 1
  }

  # Initializing the classification result
  # Note: number of columns equals the number of dates;
  # Each column refers to one date

  classification <- matrix(0, nrow = num_rows, ncol = n_dates)

  valorm <- rep(-100, num_rows)

  # Facilitating the shifting of columns in the dataframe
  index_shifts <- cumsum(c(0, n_classes[-length(n_classes)]))

  n_trajectories <- nrow(traj_weights)

  #traj_weights$final_weights <- traj_weights$pesoF

  for (s in 1:n_trajectories) {
    if (traj_weights$final_weights[s] != 0) {
      trajectory_sequence <- as.numeric(traj_weights[s, 1:n_dates])

      value <- rep(1, num_rows)

      # Multiplying the likelihood values of each class w
      # of the trajectory s
      for (w in seq_along(trajectory_sequence)) {
        col_index <- trajectory_sequence[w] + index_shifts[w]
        value <- value * data[, col_index]       # Value is now the product of the probabilities of the trajectory classes s

      }

      # Multiplying the result by the trajectory weight s
      # If using discriminant function, add...value would start with 0
      value <- value * traj_weights$final_weights[s]

      # Checking which indices have higher values than the previous one
      update_indices <- which(value > valorm)

      # Only replaces at the index where the value is greater
      # Replaces the trajectory by complet
      if (length(update_indices) > 0) {
        classification[update_indices, ] <- matrix(trajectory_sequence, nrow = length(update_indices), ncol = n_dates, byrow = TRUE)

        # Update the minimum value
        valorm[update_indices] <- value[update_indices]
      }
    }
    cat("trajectory", s,"\n")
  }

  # Writing the raster
  b1 <- setValues(raster(class_likelihood_paths[1]), rep(0, num_rows))

  for (i in 1:n_dates){

    b1[]<-classification[,i]

    if (!is.null(masking)){

      cloud_mask <- brick(masking[i])

      b1 <- masking_replace(b1,cloud_mask,0)
    }

    writeRaster(b1,paste0(output_path,"CMAP_",i,".tif"),overwrite=TRUE)

  }
}





#' Verify if All Images Have the Same Dimensions
#'
#' This function checks if all images in the provided file paths have the same dimensions.
#'
#' @param paths A character vector containing the file paths of the images to be checked.
#' @return A logical value indicating whether all images have the same dimensions (TRUE) or not (FALSE).
#' @importFrom raster raster
#' @export
verify_dimensions <- function(paths) {
  # Inicializar variáveis para armazenar as dimensões da primeira imagem
  inicial <- NULL

  # Loop para ler cada imagem e verificar suas dimensões
  for (path in paths) {
    # Carregar a imagem
    img <- raster(path)

    # Obter dimensões da imagem
    dimensoes <- c(nrow(img), ncol(img))

    # Se é a primeira imagem, definir como inicial
    if (is.null(inicial)) {
      inicial <- dimensoes
    } else {
      # Se as dimensões não forem iguais às da primeira imagem, retornar FALSE
      if (!all(inicial == dimensoes)) {
        cat("Imagem em", path, "tem dimensões diferentes: esperado", inicial, "obtido", dimensoes, "\n")
        return(FALSE)
      }
    }
  }

  # Se todas as imagens têm as mesmas dimensões
  cat("Todas as imagens têm as mesmas dimensões:", inicial, "\n")
  return(TRUE)
}

#Sintaxe-----

# output_path_C3 <- "G:/Meu Drive/INPE/projeto_dissertacao/7_CMAP/C3/"

# Importing data
# class_likelihood_paths_C3 <-c("G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/likelihoods_2005.TIF",
#                               "G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/likelihoods_2006.TIF",
#                               "G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/likelihoods_2007.TIF",
#                               "G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/likelihoods_2008.TIF",
#                               "G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/likelihoods_2009.TIF",
#                               "G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/likelihoods_2010.TIF")
#
# cloud_mask_paths_C3 <- c('G:/Meu Drive/INPE/projeto_dissertacao/0_images/3_cloud_mask/rec_cloud_mask_2005.TIF',
#                          'G:/Meu Drive/INPE/projeto_dissertacao/0_images/3_cloud_mask/rec_cloud_mask_2006.TIF',
#                          'G:/Meu Drive/INPE/projeto_dissertacao/0_images/3_cloud_mask/rec_cloud_mask_2007.TIF',
#                          'G:/Meu Drive/INPE/projeto_dissertacao/0_images/3_cloud_mask/rec_cloud_mask_2008.TIF',
#                          'G:/Meu Drive/INPE/projeto_dissertacao/0_images/3_cloud_mask/rec_cloud_mask_2009.TIF',
#                          'G:/Meu Drive/INPE/projeto_dissertacao/0_images/3_cloud_mask/rec_cloud_mask_2010.TIF')
#
# validade_C3 <- read.table("G:/Meu Drive/INPE/projeto_dissertacao/7_CMAP/C3/pesos_C3.txt",sep=";",h=T)
# peso_C3 <- validade_C3[,-c(7:12)]

# trajectory_matrix_C3 <- validade_C3

# CMAP_classifier(class_likelihood_paths_C3,trajectory_matrix_C3,output_path_C3,masking = cloud_mask_paths_C3)




# Masking


#' Replace Values in an Image Based on a Mask
#'
#' This function replaces values in an image based on a mask, applying the specified value where the mask is TRUE.
#'
#' @param image A RasterStack or RasterBrick object representing the image.
#' @param mask A RasterLayer object representing the mask.
#' @param value The value to replace in the image where the mask is TRUE.
#' @return A RasterStack or RasterBrick object with the modified values.
#' @import terra
#' @export
#'
#' @examples
#' # Carregando o pacote necessário
#' library(raster)
#'
#' # Create an example image and a mask
#' image <- stack(raster(matrix(1:100, 10, 10)), raster(matrix(1:100, 10, 10)))
#' mask <- raster(matrix(sample(0:1, 100, replace = TRUE), 10, 10))
#'
#' # Replace values based on the mask
#' result <- masking_replace(image, mask, value = 0)
#'
#' # Plot the result
#' plot(result)
masking_replace <- function(image, mask, value) {
  # Verificar se a imagem e a máscara têm o mesmo número de linhas e colunas
  if (dim(image)[1] != dim(mask)[1] || dim(image)[2] != dim(mask)[2]) {
    stop("A imagem e a máscara devem ter a mesma dimensão")
  }

  # Aplicar a máscara em todas as bandas da imagem
  for (band in 1:nlayers(image)) {
    # Obter a banda atual
    actual_band <- raster::subset(image, band)

    # Substituir os valores onde a máscara é TRUE (presumindo que TRUE indica nuvens)
    pixels <- values(actual_band)
    pixels[values(mask)==1] <- value
    values(actual_band) <- pixels

    # Salvar a banda modificada de volta na imagem
    image[[band]] <- actual_band
  }

  # Retorna a imagem modificada
  return(image)
}


#--------------------------------ASSESSMENT-------------------------------------




#' Calculate Impossible Transitions in a Time Series of Raster Classifications
#'
#' This function calculates the percentage of impossible transitions in a time series of raster classifications based on transition matrices.
#'
#' @param diretorioBase The base directory where the raster files are located.
#' @param anos A vector of years corresponding to the raster files.
#' @param caminhosTransicoes A vector of file paths to the transition matrices.
#' @return The percentage of impossible transitions.
#' @import raster
#' @export
#' @examples
#' # Assuming example classification rasters are saved in "Classification_2000.TIF", "Classification_2001.TIF"
#' diretorioBase <- "path/to/directory/"
#' anos <- c(2000, 2001)
#' caminhosTransicoes <- c("path/to/transition1.txt", "path/to/transition2.txt")
#' if (all(file.exists(paste0(diretorioBase, "Classification_", anos, ".TIF"))) && all(file.exists(caminhosTransicoes))) {
#'   porcentagemImpossiveis <- calcularTransicoesImpossiveis(diretorioBase, anos, caminhosTransicoes)
#' } else {
#'   message("Example files do not exist. Please provide the correct paths to the classification rasters and transition matrices.")
#' }
calcularTransicoesImpossiveis <- function(diretorioBase, anos, caminhosTransicoes) {
  # Inicializa a lista de rasters combinados com ponderação exponencial
  rasterCombinado <- NULL
  numTransicoes <- length(caminhosTransicoes)

  # Processar cada raster de entrada, aplicando uma ponderação exponencial decrescente
  for (i in 1:(numTransicoes + 1)) {
    # Carregar o raster correspondente ao ano
    rasterAnual <- raster(paste0(diretorioBase, "Classification_", anos[i], ".TIF"))
    # Aplicar ponderação decrescente
    rasterPonderado <- 10^(6 - i) * rasterAnual

    # Combinar rasters ponderados
    if (is.null(rasterCombinado)) {
      rasterCombinado <- rasterPonderado
    } else {
      rasterCombinado <- rasterCombinado + rasterPonderado
    }
  }

  # Armazenar informações sobre as classes de cada matriz de transição
  numClasses <- vector(length = (numTransicoes + 1))
  classes <- list(length = numTransicoes + 1)

  # Ler e armazenar dados das matrizes de transição
  for (z in 1:numTransicoes) {
    matrizTransicao <- read.table(caminhosTransicoes[z], sep = ";", header = FALSE)

    if (z == 1) {
      numClasses[1] <- dim(matrizTransicao)[1]
      classes[[1]] <- 1:dim(matrizTransicao)[1]
    }
    numClasses[z + 1] <- dim(matrizTransicao)[2]
    classes[[z + 1]] <- 1:dim(matrizTransicao)[2]
  }

  # Criar tabela de todas as combinações possíveis de transições
  tabelaTransicoes <- expand.grid(classes[])
  pesos <- vector(length = dim(tabelaTransicoes)[1])

  # Atribuir pesos às transições conforme definido nas matrizes
  for (z in 1:numTransicoes) {
    matrizTransicao <- read.table(caminhosTransicoes[z], sep = ";", header = FALSE)
    for (linha in 1:dim(matrizTransicao)[1]) {
      for (coluna in 1:dim(matrizTransicao)[2]) {
        indices <- which(tabelaTransicoes[, z] == linha & tabelaTransicoes[, z + 1] == coluna)
        pesos[indices] <- matrizTransicao[linha, coluna]
      }
    }
    tabelaTransicoes <- cbind(tabelaTransicoes, pesos)
  }

  # Calcular o produto final dos pesos para cada combinação de transição
  if (numTransicoes > 1) {
    pesosFinais <- apply(tabelaTransicoes[, ((numTransicoes + 2):dim(tabelaTransicoes)[2])], 1, prod)
    tabelaTransicoes <- cbind(tabelaTransicoes, pesoFinal = pesosFinais)
  }

  # Filtrar transições com peso final zero (impossíveis)
  transicoesImpossiveis <- tabelaTransicoes[which(tabelaTransicoes$pesoFinal == 0), 1:(numTransicoes + 1)]
  transicoesImpossiveis$codificacaoTransicoes <- 0

  # Codificar transições impossíveis para identificação única
  for (i in 1:(numTransicoes + 1)) {
    transicoesImpossiveis$codificacaoTransicoes <- transicoesImpossiveis$codificacaoTransicoes + 10^(6 - i) * transicoesImpossiveis[, i]
  }

  # Criar um raster de saída identificando transições impossíveis
  rasterFinal <- rasterCombinado * 0
  indicesImpossiveis <- which(rasterCombinado[] %in% as.vector(transicoesImpossiveis$codificacaoTransicoes))
  rasterFinal[indicesImpossiveis] <- 1

  # Escrever o raster resultante para o arquivo
  writeRaster(rasterFinal, file = paste0(diretorioBase, "trajetorias_impossiveis.tif"),overwrite = TRUE)

  # Calcular a porcentagem de transições impossíveis
  contagemTransicoes <- table(as.data.frame(rasterFinal))
  porcentagemImpossiveis <- contagemTransicoes[2] / sum(contagemTransicoes) * 100

  return(porcentagemImpossiveis)
}

# ------------C3
# caminhosTransicoes_C3  <-c("G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2005-2006.txt",
#                            "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2006-2007.txt",
#                            "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2007-2008.txt",
#                            "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2008-2009.txt",
#                            "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2009-2010.txt")
#
# diretorioBase_C3 <- "G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/"
# anos_C3 <- 2005:2010
#
# resultado_C3 <- calcularTransicoesImpossiveis(diretorioBase_C3, anos_C3, caminhosTransicoes_C3)
# print(resultado_C3) # 57.47937

# ------------C6
# caminhosTransicoes_C6  <- c("G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2005-2006.txt",
#                             "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2006-2007.txt",
#                             "G:/Meu Drive/INPE/projeto_dissertacao/2_documents/documentos_reis/CMAP1/0_matrizes/01/2008-2009.txt")
#
#
# diretorioBase_C6 <- "G:/Meu Drive/INPE/projeto_dissertacao/6_cenarios/resultados_C3/" #tenho que pegar o maxver
# anos_C6 <- c(2005,2006,2007,2009)
# resultado_C6 <- calcularTransicoesImpossiveis(diretorioBase_C6, anos_C6, caminhosTransicoes_C6)
# print(resultado_C6) # 52.20572

#' Compare Raster Classifications
#'
#' This function compares raster classifications from two sources over multiple years and calculates the differences.
#'
#' @param baseDirectory_ML The base directory for ML classification rasters.
#' @param baseDirectory_CMAP The base directory for CMAP classification rasters.
#' @param years A vector of years corresponding to the raster files.
#' @param maskPath The path to the cloud mask raster file.
#' @return A list containing the accumulated differences and the final percentage of disagreement.
#' @import raster
#' @export
#' @examples
#' # Assuming example classification rasters are saved in "ML/Classification_2000.TIF", "CMAP/CMAP_1.TIF"
#' baseDirectory_ML <- "path/to/ML/"
#' baseDirectory_CMAP <- "path/to/CMAP/"
#' years <- c(2000, 2001, 2002)
#' maskPath <- "path/to/cloud_mask.TIF"
#' if (file.exists(maskPath) && all(file.exists(paste0(baseDirectory_ML, "Classification_", years, ".TIF"))) && all(file.exists(paste0(baseDirectory_CMAP, "CMAP_", 1:length(years), ".TIF")))) {
#'   result <- compareRasterClassification(baseDirectory_ML, baseDirectory_CMAP, years, maskPath)
#' } else {
#'   message("Example files do not exist. Please provide the correct paths to the classification rasters and cloud mask.")
#' }
compareRasterClassification <- function(baseDirectory_ML, baseDirectory_CMAP, years, maskPath) {

  cloudMask <- raster(maskPath)
  # Initialization of a variable to accumulate differences across years
  accumulatedDifferences <- NULL

  # Loop over each year to process the corresponding rasters
  for (i in 1:length(years)) {
    # Construct paths to the ML and CMAP rasters
    mlRasterPath <- paste0(baseDirectory_ML,"Classification_", years[i], ".TIF")
    cmapRasterPath <- paste0(baseDirectory_CMAP,  "CMAP_", i, ".TIF")

    # Load the rasters
    mlRaster <- raster(mlRasterPath)
    cmapRaster <- raster(cmapRasterPath)

    # Calculate differences
    difference <- mlRaster - cmapRaster
    differenceMask <- mlRaster * 0
    differenceMask[which(difference[] != 0)] <- 1

    # Apply cloud mask to the difference mask
    cloudAdjustedDifference <- differenceMask
    cloudAdjustedDifference[which(cloudMask[] == 1)] <- 2

    # Calculate and store the difference tables for each year
    differenceTable <- table(as.data.frame(difference))
    adjustedDifferenceTable <- table(as.data.frame(cloudAdjustedDifference))

    # Accumulate the differences in a trajectory raster
    if (is.null(accumulatedDifferences)) {
      accumulatedDifferences <- differenceMask
    } else {
      accumulatedDifferences <- accumulatedDifferences + differenceMask
    }

    # Print current year and percentage of significant differences
    print(years[i])
    if ("1" %in% rownames(adjustedDifferenceTable)) {
      print(100 * adjustedDifferenceTable["1"] / (adjustedDifferenceTable["0"] + adjustedDifferenceTable["1"]))
    } else {
      print(0)
    }
  }

  # Adjust accumulated differences with cloud mask
  accumulatedDifferences[which(cloudMask[] == 1)] <- -1
  finalDifferenceTable <- table(as.data.frame(accumulatedDifferences))

  # Print final trajectory analysis
  cat("Total disagreement:")
  if ("1" %in% rownames(finalDifferenceTable)) {
    print(100 * sum(finalDifferenceTable[rownames(finalDifferenceTable) > 0]) / sum(finalDifferenceTable[rownames(finalDifferenceTable) > -1]))
  } else {
    print(0)
  }

  return(list(accumulatedDifferences = accumulatedDifferences, finalPercentage = finalDifferenceTable))
}
