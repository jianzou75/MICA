library(parallel)

Simulation <- function(setting){
  study.data.generator <- function(mu.matrix, sd.vector, grp.size){
    study_data  <- rbind(
      cbind(matrix(rnorm(grp.size[1]*300,  mu.matrix[1,1], sd.vector[1]), nrow = 300),  matrix(rnorm(grp.size[2]*300,  mu.matrix[1,2], sd.vector[1]), nrow = 300),  matrix(rnorm(grp.size[3]*300,  mu.matrix[1,3], sd.vector[1]), nrow = 300)), ## cateogry 1
      cbind(matrix(rnorm(grp.size[1]*100,  mu.matrix[2,1], sd.vector[2]), nrow = 100),  matrix(rnorm(grp.size[2]*100,  mu.matrix[2,2], sd.vector[2]), nrow = 100),  matrix(rnorm(grp.size[3]*100,  mu.matrix[2,3], sd.vector[2]), nrow = 100)), ## cateogry 2
      cbind(matrix(rnorm(grp.size[1]*100,  mu.matrix[3,1], sd.vector[3]), nrow = 100),  matrix(rnorm(grp.size[2]*100,  mu.matrix[3,2], sd.vector[3]), nrow = 100),  matrix(rnorm(grp.size[3]*100,  mu.matrix[3,3], sd.vector[3]), nrow = 100)), ## cateogry 3
      cbind(matrix(rnorm(grp.size[1]*1500, mu.matrix[4,1], sd.vector[4]), nrow = 1500), matrix(rnorm(grp.size[2]*1500, mu.matrix[4,2], sd.vector[4]), nrow = 1500), matrix(rnorm(grp.size[3]*1500, mu.matrix[4,3], sd.vector[4]), nrow = 1500)) ## cateogry 4
    )
    study_label <- c(rep(1, grp.size[1]), rep(2, grp.size[2]), rep(3, grp.size[3]))
    return(list(study_data = study_data,
                study_label = study_label))
  }


  if(setting == "case1"){ ## Effect size = 0.5
    study1.mu  <- matrix(c(1,3,5, 1,3,5, 1,3,5, 0,0,0), nrow = 4, byrow = T)
    study1.sd  <- rep(3.5, 4)
    study1.grp <- c(10, 5, 8)
    study1 <- study.data.generator(mu.matrix = study1.mu, sd.vector = study1.sd, grp.size = study1.grp)

    study2.mu  <- matrix(c(2,4,6, 6,4,2, 2,4,6, 0,0,0), nrow = 4, byrow = T)
    study2.sd  <- rep(3.1, 4)
    study2.grp <- c(5,8,10)
    study2 <- study.data.generator(mu.matrix = study2.mu, sd.vector = study2.sd, grp.size = study2.grp)

    study3.mu  <- matrix(c(1,4,7, 1,7,1, 0,0,0, 0,0,0), nrow = 4, byrow = T)
    study3.sd  <- c(4.4, 5.9, 4.4, 4.4)
    study3.grp <- c(8,10,5)
    study3 <- study.data.generator(mu.matrix = study3.mu, sd.vector = study3.sd, grp.size = study3.grp)
  }

  if(setting == "case2"){ ## Effect size = 0.5
    study1.mu  <- matrix(c(1,3,5, 1,3,5, 1,3,5, 0,0,0), nrow = 4, byrow = T)
    study1.sd  <- rep(2.9, 4)
    study1.grp <- c(10, 5, 8)
    study1 <- study.data.generator(mu.matrix = study1.mu, sd.vector = study1.sd, grp.size = study1.grp)

    study2.mu  <- matrix(c(2,4,6, 6,4,2, 2,4,6, 0,0,0), nrow = 4, byrow = T)
    study2.sd  <- rep(2.6, 4)
    study2.grp <- c(5,8,10)
    study2 <- study.data.generator(mu.matrix = study2.mu, sd.vector = study2.sd, grp.size = study2.grp)

    study3.mu  <- matrix(c(1,4,7, 1,7,1, 0,0,0, 0,0,0), nrow = 4, byrow = T)
    study3.sd  <- c(3.7, 4.8, 3.7, 3.7)
    study3.grp <- c(8,10,5)
    study3 <- study.data.generator(mu.matrix = study3.mu, sd.vector = study3.sd, grp.size = study3.grp)
  }

  if(setting == "case3"){ ## Effect size = 0.5
    study1.mu  <- matrix(c(1,3,5, 1,3,5, 1,3,5, 0,0,0), nrow = 4, byrow = T)
    study1.sd  <- rep(2.5, 4)
    study1.grp <- c(10, 5, 8)
    study1 <- study.data.generator(mu.matrix = study1.mu, sd.vector = study1.sd, grp.size = study1.grp)

    study2.mu  <- matrix(c(2,4,6, 6,4,2, 2,4,6, 0,0,0), nrow = 4, byrow = T)
    study2.sd  <- rep(2.2, 4)
    study2.grp <- c(5,8,10)
    study2 <- study.data.generator(mu.matrix = study2.mu, sd.vector = study2.sd, grp.size = study2.grp)

    study3.mu  <- matrix(c(1,4,7, 1,7,1, 0,0,0, 0,0,0), nrow = 4, byrow = T)
    study3.sd  <- c(3.2, 4.3, 3.2, 3.2)
    study3.grp <- c(8,10,5)
    study3 <- study.data.generator(mu.matrix = study3.mu, sd.vector = study3.sd, grp.size = study3.grp)
  }

  return(list(study.data.list = c(list(t(study1$study_data)), list(t(study2$study_data)), list(t(study3$study_data))),
              study.label.list = c(list(study1$study_label), list(study2$study_label), list(study3$study_label))))

}


## Computation
### Case 1
perm.test.res.list1 <- list()
for(i in 1:200){
  set.seed(i)
  sim.data = Simulation("case1")
  perm.test.res <- mc.permute.test(study.data.list = sim.data$study.data.list, study.label.list = sim.data$study.label.list,
                                   w.est = NULL, n.perm = 1000, n.parallel = 50)
  perm.test.res.list1 <- c(perm.test.res.list1, list(perm.test.res))
}



### Case 2
perm.test.res.list2 <- list()
for(i in 1:200){
  set.seed(i)
  sim.data = Simulation("case2")
  perm.test.res <- mc.permute.test(study.data.list = sim.data$study.data.list, study.label.list = sim.data$study.label.list,
                                   w.est = NULL, n.perm = 1000, n.parallel = 50)
  perm.test.res.list2 <- c(perm.test.res.list2, list(perm.test.res))
}

### Case 3
perm.test.res.list3 <- list()
for(i in 1:200){
  set.seed(i)
  sim.data = Simulation("case3")
  perm.test.res <- mc.permute.test(study.data.list = sim.data$study.data.list, study.label.list = sim.data$study.label.list,
                                   w.est = NULL, n.perm = 1000, n.parallel = 50)
  perm.test.res.list3 <- c(perm.test.res.list3, list(perm.test.res))
}



### After iteration
mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$TC_tbl$qval <= 0.05)[1:300])))
mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$minMI_tbl$qval <= 0.05)[1:300])))

mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$TC_tbl$qval <= 0.05)[301:400])))
mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$minMI_tbl$qval <= 0.05)[301:400])))

mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$TC_tbl$qval <= 0.05)[401:500])))
mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$minMI_tbl$qval <= 0.05)[401:500])))

mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$TC_tbl$qval <= 0.05)[501:2000])))
mean(sapply(1:50, function(x) sum((perm.test.res.list[[x]]$minMI_tbl$qval <= 0.05)[501:2000])))






