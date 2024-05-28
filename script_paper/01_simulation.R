## Simulation setting
Simulation <- function(setting){
  study.data.generator <- function(mu.matrix, sd.vector, grp.size){
    study_data  <- c(rnorm(grp.size[1], mu.matrix[1], sd.vector[1]),
                     rnorm(grp.size[2], mu.matrix[2], sd.vector[2]),
                     rnorm(grp.size[3], mu.matrix[3], sd.vector[3]))
    study_label <- c(rep(1, grp.size[1]), rep(2, grp.size[2]), rep(3, grp.size[3]))
    return(list(study_data = study_data,
                study_label = study_label))
  }


  if(setting == "case1"){ ## complelety concordant
    study1.mu  <- c(1, 3, 5)
    study1.sd  <- rep(5, 3)
    study1.grp <- c(10, 10, 10)
    study1 <- study.data.generator(mu.matrix = study1.mu, sd.vector = study1.sd, grp.size = study1.grp)

    study2.mu  <- c(1, 3, 5)
    study2.sd  <- rep(5, 3)
    study2.grp <- c(10, 10, 10)
    study2 <- study.data.generator(mu.matrix = study2.mu, sd.vector = study2.sd, grp.size = study2.grp)

    study3.mu  <- c(1, 3, 5)
    study3.sd  <- rep(5, 3)
    study3.grp <- c(10, 10, 10)
    study3 <- study.data.generator(mu.matrix = study3.mu, sd.vector = study3.sd, grp.size = study3.grp)

    study4.mu  <- c(1, 3, 5)
    study4.sd  <- rep(5, 3)
    study4.grp <- c(10, 10, 10)
    study4 <- study.data.generator(mu.matrix = study4.mu, sd.vector = study4.sd, grp.size = study4.grp)
  }

  if(setting == "case2"){ ## Effect size = 0.5
    study1.mu  <- c(5, 3, 1)
    study1.sd  <- rep(5, 3)
    study1.grp <- c(10, 10, 10)
    study1 <- study.data.generator(mu.matrix = study1.mu, sd.vector = study1.sd, grp.size = study1.grp)

    study2.mu  <- c(5, 3, 1)
    study2.sd  <- rep(5, 3)
    study2.grp <- c(10, 10, 10)
    study2 <- study.data.generator(mu.matrix = study2.mu, sd.vector = study2.sd, grp.size = study2.grp)

    study3.mu  <- c(5, 3, 1)
    study3.sd  <- rep(5, 3)
    study3.grp <- c(10, 10, 10)
    study3 <- study.data.generator(mu.matrix = study3.mu, sd.vector = study3.sd, grp.size = study3.grp)

    study4.mu  <- c(1, 7, 1)
    study4.sd  <- rep(5, 3)
    study4.grp <- c(10, 10, 10)
    study4 <- study.data.generator(mu.matrix = study4.mu, sd.vector = study4.sd, grp.size = study4.grp)
  }

  if(setting == "case3"){ ## Effect size = 0.5
    study1.mu  <- c(1, 3, 5)
    study1.sd  <- rep(5, 3)
    study1.grp <- c(10, 10, 10)
    study1 <- study.data.generator(mu.matrix = study1.mu, sd.vector = study1.sd, grp.size = study1.grp)

    study2.mu  <- c(1, 3, 5)
    study2.sd  <- rep(5, 3)
    study2.grp <- c(10, 10, 10)
    study2 <- study.data.generator(mu.matrix = study2.mu, sd.vector = study2.sd, grp.size = study2.grp)

    study3.mu  <- c(1, 7, 1)
    study3.sd  <- rep(5, 3)
    study3.grp <- c(10, 10, 10)
    study3 <- study.data.generator(mu.matrix = study3.mu, sd.vector = study3.sd, grp.size = study3.grp)

    study4.mu  <- c(1, 7, 1)
    study4.sd  <- rep(5, 3)
    study4.grp <- c(10, 10, 10)
    study4 <- study.data.generator(mu.matrix = study4.mu, sd.vector = study4.sd, grp.size = study4.grp)
  }

  if(setting == "case4"){ ## Effect size = 0.5
    study1.mu  <- c(0, 0, 0)
    study1.sd  <- rep(5, 3)
    study1.grp <- c(10, 10, 10)
    study1 <- study.data.generator(mu.matrix = study1.mu, sd.vector = study1.sd, grp.size = study1.grp)

    study2.mu  <- c(0, 0, 0)
    study2.sd  <- rep(5, 3)
    study2.grp <- c(10, 10, 10)
    study2 <- study.data.generator(mu.matrix = study2.mu, sd.vector = study2.sd, grp.size = study2.grp)

    study3.mu  <- c(0, 0, 0)
    study3.sd  <- rep(5, 3)
    study3.grp <- c(10, 10, 10)
    study3 <- study.data.generator(mu.matrix = study3.mu, sd.vector = study3.sd, grp.size = study3.grp)

    study4.mu  <- c(0, 0, 0)
    study4.sd  <- rep(5, 3)
    study4.grp <- c(10, 10, 10)
    study4 <- study.data.generator(mu.matrix = study4.mu, sd.vector = study4.sd, grp.size = study4.grp)
  }

  return(list(study.data.list = c(list(study1$study_data), list(study2$study_data), list(study3$study_data), list(study4$study_data)),
              study.label.list = c(list(study1$study_label), list(study2$study_label), list(study3$study_label), list(study4$study_label))))

}

## Module 1
module1.data <- list()
module1.result <- list()
for(i in 1:500){
  set.seed(i)
  sim.data <- Simulation("case1")
  res <- mica(sim.data$study.data.list, sim.data$study.label.list, post.hoc = FALSE)
  module1.data <- append(module1.data, list(sim.data))
  module1.result <- append(module1.result, list(res))
}

## Module 2
module2.data <- list()
module2.result <- list()
for(i in 1:500){
  set.seed(i)
  sim.data <- Simulation("case2")
  res <- mica(sim.data$study.data.list, sim.data$study.label.list, post.hoc = FALSE)
  module2.data <- append(module2.data, list(sim.data))
  module2.result <- append(module2.result, list(res))
}

## Module 3
module3.data <- list()
module3.result <- list()
for(i in 1:500){
  set.seed(i)
  sim.data <- Simulation("case3")
  res <- mica(sim.data$study.data.list, sim.data$study.label.list, post.hoc = FALSE)
  module3.data <- append(module3.data, list(sim.data))
  module3.result <- append(module3.result, list(res))
}

## Module 4
module4.data <- list()
module4.result <- list()
for(i in 1:500){
  set.seed(i)
  sim.data <- Simulation("case4")
  res <- mica(sim.data$study.data.list, sim.data$study.label.list, post.hoc = FALSE)
  module4.data <- append(module4.data, list(sim.data))
  module4.result <- append(module4.result, list(res))
}

## Genome wide analysis
genome_matrix_simu <- function(case){
  data_list <- lapply(1:500, function(x) Simulation(case)$study.data.list)
  study1 <- t(sapply(1:500, function(x) data_list[[x]][[1]]))
  study2 <- t(sapply(1:500, function(x) data_list[[x]][[2]]))
  study3 <- t(sapply(1:500, function(x) data_list[[x]][[3]]))
  study4 <- t(sapply(1:500, function(x) data_list[[x]][[4]]))
  return(list(study1 = study1, study2 = study2, study3 = study3, study4 = study4))
}
for(i in 1:200){
  set.seed(i)
  module1 <- genome_matrix_simu("case1")
  module2 <- genome_matrix_simu("case2")
  module3 <- genome_matrix_simu("case3")
  module4 <- genome_matrix_simu("case4")
  study1 <- rbind(module1$study1, module2$study1, module3$study1, module4$study1)
  study2 <- rbind(module1$study2, module2$study2, module3$study2, module4$study2)
  study3 <- rbind(module1$study3, module2$study3, module3$study3, module4$study3)
  study4 <- rbind(module1$study4, module2$study4, module3$study4, module4$study4)
  sim.data <- c(list(study1), list(study2), list(study3), list(study4))
  sim.label <- Simulation("case1")$study.label.list
  sim_result <- mica.full(sim.data, sim.label, post.hoc = FALSE)
  genome.result <- append(genome.result, list(sim_result[1:4]))
}



