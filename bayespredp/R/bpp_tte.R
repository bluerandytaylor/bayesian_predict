library(survival)
library(R2jags)
library(snowfall)

# 重要函数 ---------------------------------------------------------------------
simes.test <- function(x){
  r=rank(x,ties.method = "random")
  T=min(length(x)*x/r)
  return(T)
}

# 1.场景设置 -------------------------------------------------------------------

## 1.1.疗效参数 ----------------------------------------------------------------
#考虑有3个试验臂，一个对照臂。
sc_binary.vec <- c(0.2,0.2,0.2,0.2,#1
                   0.2,0.3,0.35,0.4,#2
                   0.2,0.2,0.2,0.2,#3
                   0.2,0.2,0.2,0.2,#4
                   0.2,0.3,0.35,0.4,#5
                   0.2,0.3,0.375,0.425,#6
                   0.2,0.2,0.2,0.2,#7
                   0.2,0.3,0.35,0.4)#8
sc_binary <-  matrix(sc_binary.vec,nrow = length(sc_binary.vec)/4,byrow=TRUE)
colnames(sc_binary) <- c("pc","p1","p2","p3")

siglevel <- 0.025

#基线风险服从weibull分布，尺度参数为lambda，形状参数为gamma。lambda直接乘在最前面。
#notice:需要进行调整，使得期中分析时有合适的右删失比例
#与working plan中的基线风险表达不同
#这两个字母来源，包括下面的survival time的生成，应该是来自于：
#https://cran.r-project.org/web/packages/rsimsum/vignettes/B-relhaz.html
nu_wp <- 2
lambda_wp <- 7

gamma <- nu_wp
lambda <- 1/((lambda_wp)^nu_wp)

sc_hazard.vec <- c(0,0,0,#1
                   0,0,0,#2
                   0,0,0,#3
                   -0.2,-0.35,-0.45,#4
                   -0.2,-0.35,-0.45,#5
                   0,0,0,#6
                   -0.2,-0.35,-0.45,#7
                   -0.1,-0.2,-0.25)#8
sc_hazard <- matrix(sc_hazard.vec,nrow = length(sc_hazard.vec)/3,byrow = TRUE)
colnames(sc_hazard) <- c("beta1","beta2","beta3")
#beta1,2,3分别为三个试验药物对疗效的直接作用

sc_indicate <- c(0,0,-1,0,0,-2,-1,-1)#为zeta的取值，代表8个场景下的短期终点对长期终点的预测作用

sc_hr1 <- exp(sc_hazard[,"beta1"]+(sc_binary[,"p1"]-sc_binary[,"pc"])*sc_indicate)
sc_hr2 <- exp(sc_hazard[,"beta2"]+(sc_binary[,"p2"]-sc_binary[,"pc"])*sc_indicate)
sc_hr3 <- exp(sc_hazard[,"beta3"]+(sc_binary[,"p3"]-sc_binary[,"pc"])*sc_indicate)
hrs <- data.frame(hr1=sc_hr1,
                  hr2=sc_hr2,
                  hr3=sc_hr3)
## 1.2. 时间点与两阶段样本量 ---------------------------------------------------
tau <- 13#期中分析时间点
L <- 24#患者招募最长时间
t_end <- 36#试验结束时间

ss_eacharm <- 100#暂定为既定的整个试验的样本量。

#入组时间。需要注意ss_eacharm需要能够被25整除，否则会报错。
t_join <- rep(rep(seq(0,24,length=25),each=ss_eacharm/25),ncol(sc_binary))
#第一阶段每个臂上患者数量
ss_st1 <- sum(t_join<tau)/ncol(sc_binary)
#预计第二阶段每个臂上患者数量
ss_st2_plan <- ss_eacharm-ss_st1
#权重计算
scontrol <- function(x){exp(-(x/lambda_wp)^nu_wp)}
st <- function(x){exp(-((x/lambda_wp)^nu_wp)*0.6)}
peventsc <- 1-(scontrol(t_end-tau+1)+4*scontrol(t_end-tau/2+1/2)+scontrol(t_end))/6
peventst <- 1-(st(t_end-tau+1)+4*st(t_end-tau/2+1/2)+st(t_end))/6

pmean <- mean(c(peventsc,peventst))

omega1 <- ((ss_st1*pmean)/ss_eacharm)^(1/2)
omega2 <- (1-(ss_st1*pmean)/ss_eacharm)^(1/2)

chw <- function(s1p,s2p){
  pchw <- 1-pnorm(omega1*qnorm(1-s1p)+omega2*qnorm(1-s2p))
  return(pchw)
}

#样本量再估计时需考虑，两个阶段最大样本量
ss_eacharm_max <- 150
#样本量再估计时需考虑，第二阶段最大样本量
ss_st2_max <- ss_eacharm_max-ss_st1

## 1.3. 期中决策参数 -----------------------------------------------------------

#favorable zone
epsilon <- 0.1
#promising zone
phi <- 0.6
#unfavorable zone
eta <- 0.1
#futility zone


  ## 2.1.患者的治疗臂数据 --------------------------------------------------------
  
  trt <- rbind(matrix(c(rep(0,ss_eacharm*3)),byrow = FALSE,ncol = 3),#对照臂
               matrix(c(rep(1,ss_eacharm),rep(0,ss_eacharm*2)),byrow = FALSE,ncol = 3),#试验臂1
               matrix(c(rep(0,ss_eacharm),rep(1,ss_eacharm),rep(0,ss_eacharm)),byrow=FALSE,ncol=3),#试验臂2
               matrix(c(rep(0,ss_eacharm*2),rep(1,ss_eacharm)),byrow = FALSE,ncol = 3))#试验臂3
  
  ## 2.2.二分类终点数据 -----------------------------------------------------------
  
  binary <- c(rbinom(ss_eacharm,1,rates[1]),#对照臂
              rbinom(ss_eacharm,1,rates[2]),#试验臂1
              rbinom(ss_eacharm,1,rates[3]),#试验臂2
              rbinom(ss_eacharm,1,rates[4]))#试验臂3
  
  ## 2.3.长期终点数据 -------------------------------------------------------------
  
  #https://cran.r-project.org/web/packages/rsimsum/vignettes/B-relhaz.html
  u <- runif(4*ss_eacharm)
  t_surv <- abs((-log(u)/(lambda*exp(trt%*%as.matrix(hazard.pars)+binary*pred)))^(1/gamma))#寿命长度
  t_death <- t_join+t_surv#真实的死亡时间
  
  ## 2.4.期中分析整理数据 ---------------------------------------------------------
  
  t_interim_obs <- pmin(t_death,tau)
  interim_censor <- ifelse(t_death<=tau,0,1)#0表示期中分析时未删失，1表示删失
  
  #合并数据，生成原始数据
  data <- data.frame(join_time=t_join,
                     death_time=t_death,
                     interim_time=t_interim_obs-t_join,#期中分析时观察到的生存时间长度？
                     interim_obs_status=1-interim_censor,#1表示观测到死亡，status=0表示删失。
                     response=binary)
  
  #对上述数据进行筛选，只保留期中数据。
  interim_data <- data[t_join<tau,]
  interim_trt <- trt[t_join<tau,]
  interimdf <- cbind(interim_trt,interim_data)
  
  interimdf$arm[interimdf$`1`==0&interimdf$`2`==0&interimdf$`3`==0] <- "control"
  interimdf$arm[interimdf$`1`==1] <- "treatment1"
  interimdf$arm[interimdf$`2`==1] <- "treatment2"
  interimdf$arm[interimdf$`3`==1] <- "treatment3"
  
  
  ia_censorrate_arm1 <- 1-mean(interimdf$interim_obs_status[interimdf$arm=="treatment1"])
  ia_censorrate_arm2 <- 1-mean(interimdf$interim_obs_status[interimdf$arm=="treatment2"])
  ia_censorrate_arm3 <- 1-mean(interimdf$interim_obs_status[interimdf$arm=="treatment3"])
  ia_censorrate_contr <- 1-mean(interimdf$interim_obs_status[interimdf$arm=="control"])
  
  
  
  
  ## 2.5.期中分析后验计算 --------------------------------------------------------
  
  ### 2.5.1.r2jags数据 -----------------------------------------------------------
  
  n <- length(interimdf[,1])
  t <- interimdf$interim_time
  is.na(t) <- interimdf$interim_obs_status==0 
  is.censored <- 1-interimdf$interim_obs_status#是否删失，1表示删失，0表示未删，观测到死亡结局
  t.cen <- interimdf$interim_time+interimdf$interim_obs_status#若删失，则为删失时间，若未删失，则必须大于观测时间。
  trt_jags <- interim_trt
  response <- interimdf$response
  
  rjagsdata <- list(n=n,t=t,is.censored=is.censored,t.cen=t.cen,trt=trt_jags,response=response)
  
  ### 2.5.2.jags model ----------------------------------------------------------
  
  model <- "
model{
    #likehood of binary endpoint
    for (i in 1:n){
        response[i]~dbinom(p[i],1)
        p[i] <- pnorm(gamma0+gamma1*trt[i,1]+gamma2*trt[i,2]+gamma3*trt[i,3],0,1)
    }
    #likehood of tte with right censoring
    for (i in 1:n){
        is.censored[i]~dinterval(t[i],t.cen[i])
        t[i]~dweib(alpha,lambda[i])
        lambda[i] <- exp(beta0+beta1*trt[i,1]+beta2*trt[i,2]+beta3*trt[i,3]+zeta*response[i])
    }
    #priors 
    gamma0 ~ dnorm(0,0.001)
    gamma1 ~ dnorm(0,0.001)
    gamma2 ~ dnorm(0,0.001)
    gamma3 ~ dnorm(0,0.001)
  
    alpha~dgamma(1.1,1.1)
    beta0 ~ dnorm(0,0.001)
    beta1 ~ dnorm(0,0.001)
    beta2 ~ dnorm(0,0.001)
    beta3 ~ dnorm(0,0.001)
    zeta ~ dnorm(0,0.001)
}
"

parameters <- c("gamma0","gamma1","gamma2","gamma3",
                "alpha","beta0","beta1","beta2","beta3","zeta")
#模型
interimjag <- jags.model(textConnection(model),data=rjagsdata,n.chains=1)
#burn in
update(interimjag,5000)
posamples <- coda.samples(interimjag,parameters,10000)
#转换
df <- do.call(rbind.data.frame,posamples) 

# output <- jags(data = rjagsdata, model.file=textConnection(model), parameters.to.save=parameters, n.burnin = 5000,
#                n.chains = 1, n.iter = 10000, n.thin = 1,  jags.module = c("dic"),progress.bar = "none",quiet =TRUE)
# output

## 2.6.后验预测概率函数 -----------------------------------------------------


#### 2.6.1 基于单组后验，模拟是否成功 ----------------------------------------------------

single_postsamp_predict <- function(targetarm,
                                    alpha,
                                    beta0,
                                    beta1,
                                    beta2,
                                    beta3,
                                    zeta,
                                    gamma0,
                                    gamma1,
                                    gamma2,
                                    gamma3,
                                    n2plan){

  
  #提取特定试验臂与对照臂数据
  armdf1 <- interimdf[interimdf$arm %in% c("control",targetarm),]
  
  
  #### 2.6.1.1.type1:IA时已入组但删失 -------------------------------------------------
  
  type1 <- armdf1[armdf1$interim_obs_status==0,]
  #计算该类患者的生存概率上限
  type1$survuplimit <- exp(-(type1$interim_time)^alpha*exp(beta0+beta1*type1$`1`+beta2*type1$`2`+beta3*type1$`3`+zeta*type1$response))
  #从均匀分布中抽样
  type1$unifsamp <- runif(nrow(type1),min = rep(0,nrow(type1)),max=type1$survuplimit)
  #均匀分布转化为生存时间
  # type1$survt <- (-log(type1$unifsamp)/exp(beta0+beta1*type1$`1`+beta2*type1$`2`+beta3*type1$`3`+zeta*type1$response))^alpha
  #0313修改：认为计算有误
  type1$survt <- (-log(type1$unifsamp)/exp(beta0+beta1*type1$`1`+beta2*type1$`2`+beta3*type1$`3`+zeta*type1$response))^(1/alpha)
  #计算最终分析时的生存时间长度，以及是否删失
  type1$t_final <- pmin((type1$join_time+type1$survt),t_end)-type1$join_time
  type1$censor_final <- ifelse(type1$join_time+type1$survt<=t_end,0,1)#死亡为0，删失为1
  
  
  #### 2.6.1.2.type2:未入组 -------------------------------------------------------
  
  #针对选择的不同的试验臂，生成不同的治疗矩阵。前n2为对照臂，后n2为试验臂
  if (targetarm=="treatment1"){
    trt_type2 <- rbind(matrix(rep(0,n2plan*3),byrow = FALSE,ncol = 3),
                       matrix(c(rep(1,n2plan),rep(0,n2plan*2)),byrow = FALSE,ncol = 3))
  } else if(targetarm=="treatment2"){
    trt_type2 <- rbind(matrix(rep(0,n2plan*3),byrow = FALSE,ncol = 3),
                       matrix(c(rep(0,n2plan),rep(1,n2plan),rep(0,n2plan)),byrow = FALSE,ncol = 3))
  } else if(targetarm=="treatment3"){
    trt_type2 <- rbind(matrix(rep(0,n2plan*3),byrow = FALSE,ncol = 3),
                       matrix(c(rep(0,n2plan),rep(0,n2plan),rep(1,n2plan)),byrow = FALSE,ncol = 3))
  }
  #生成入组时间
  t_recruit <- runif(n2plan*2,min = rep(tau,n2plan*2),max = rep(L,n2plan*2))
  #模拟生成二分类终点,前n为对照臂，后n为试验臂
  response_type2 <- rbinom(n2plan*2,1,pnorm(gamma0+trt_type2%*%c(gamma1,gamma2,gamma3)))
  #基于二分类终点，生成生存时间，前n为对照臂，后n为试验臂(notice,r中weibull分布参数可能不同)
  # survt_type2 <- rweibull(n=n2plan*2,shape = alpha,scale = exp(-beta0-trt_type2%*%c(beta1,beta2,beta3)-zeta*response_type2))
  #notice:0312对生存时间的生成进行修改。
  u_type2 <- runif(2*n2plan)
  survt_type2 <- abs(-log(u_type2)/(exp(beta0+trt_type2%*%c(beta1,beta2,beta3)+zeta*response_type2)))^(1/alpha)
  #计算最终分析时的时间，以及是否删失。notice，这里可能有问题，type2里面删失太少了
  t_final_type2 <- pmin(survt_type2+t_recruit,t_end)-t_recruit
  censor_final_type2 <- ifelse(survt_type2+t_recruit<=t_end,0,1)#死亡为0，删失为1
  
  #### 2.6.1.3.合并所有types -------------------------------------------------------
  
  #将所有的受试者的数据全都合并起来，进行log-rank检验，判断药物是否有效。
  interim_notcensor <- interimdf[interimdf$interim_obs_status==1,]
  
  data_type0 <- data.frame(treat1=interim_notcensor$`1`,
                           treat2=interim_notcensor$`2`,
                           treat3=interim_notcensor$`3`,
                           arm=interim_notcensor$arm,
                           t_obs=interim_notcensor$interim_time,
                           censor=0)
  
  data_type1 <- data.frame(treat1=type1$`1`,
                           treat2=type1$`2`,
                           treat3=type1$`3`,
                           arm=type1$arm,
                           t_obs=type1$t_final,
                           censor=type1$censor_final)
  
  data_type2 <- data.frame(treat1=trt_type2[,1],
                           treat2=trt_type2[,2],
                           treat3=trt_type2[,3],
                           t_obs=t_final_type2,
                           censor=censor_final_type2)
  
  data_type2$arm[data_type2$treat1==0&data_type2$treat2==0&data_type2$treat3==0] <- "control"
  data_type2$arm[data_type2$treat1==1] <- "treatment1"
  data_type2$arm[data_type2$treat2==1] <- "treatment2"
  data_type2$arm[data_type2$treat3==1] <- "treatment3"
  
  #合并数据
  alldata <- rbind(data_type0,data_type1,data_type2)
  #只保留某一个试验臂与对照比数据
  alldata_singletreat <- alldata[alldata$arm %in% c("control",targetarm),]
  #status: censoring status 1=censored, 2=dead
  alldata_singletreat$status <- 2-alldata_singletreat$censor
  #https://mp.weixin.qq.com/s/19zna0b20WTudm9xgbjT0w
  
  surv_diff <- survdiff(Surv(t_obs,status)~arm,data = alldata_singletreat)
  pval <- 1-pchisq(surv_diff$chisq,length(surv_diff$n)-1)
  
  betterarm <- names(surv_diff$n)[surv_diff$obs/surv_diff$exp==min(surv_diff$obs/surv_diff$exp)]
  
  if (betterarm!="arm=control"){onesidep <- pval/2}
  if (betterarm=="arm=control"){onesidep <- 1-pval/2}
  
  return(onesidep)
  
}


### 2.6.2.mapply获得后验预测概率 ------------------------------------------------------------

bpp_predict <- function(target_arm,
                        postsamplesdf,
                        n2,
                        siglevel){
  pvals <- mapply(single_postsamp_predict,
                  postsamplesdf$alpha,
                  postsamplesdf$beta0,
                  postsamplesdf$beta1,
                  postsamplesdf$beta2,
                  postsamplesdf$beta3,
                  postsamplesdf$zeta,
                  postsamplesdf$gamma0,
                  postsamplesdf$gamma1,
                  postsamplesdf$gamma2,
                  postsamplesdf$gamma3,
                  MoreArgs = list(targetarm=target_arm,n2plan=n2)
  )
  bppvalue <- mean(pvals<siglevel)
  
  return(bppvalue)
  
}

## 2.7.SSR函数 ---------------------------------------------------------------

ssr_bpp <- function(target_arm_ssr,
                    postsamplesdf_ssr,
                    siglevel_ssr,
                    n2st_ssr,
                    n2max_ssr){
  size <- n2st_ssr
  repeat{
    bayes_pred_p <-bpp_predict(target_arm = target_arm_ssr,
                               postsamplesdf = postsamplesdf_ssr,
                               n2=size,
                               siglevel = siglevel_ssr)
    
    if (1-epsilon>bayes_pred_p&bayes_pred_p>=1-epsilon-0.1){size <- size+10}
    if (1-epsilon-0.1>bayes_pred_p&bayes_pred_p>=1-epsilon-0.2){size <- size+20}
    if (1-epsilon-0.2>bayes_pred_p){size <- size+30}
    if (bayes_pred_p>=1-epsilon){break}
    if (size>n2max_ssr){break}
    cat("size=",size," bayes_pred_p=",bayes_pred_p,"\n")
  }
  if (size>=n2max_ssr){size <- n2max_ssr}
  return(size)
}

