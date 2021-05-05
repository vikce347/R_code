myCSVfile <- "C:/Users/vikto/Desktop/Bakalaurinis/MOD.csv" #žalų trikampio nuskaitymas iš failo
    dat <- read.csv(file=myCSVfile)
    library(ChainLadder)
    library(MonteCarlo)
    library(truncnorm)
    tri <- dat[,-1]
    tri <- as.matrix(tri)
    dimnames(tri) <- list(origin=dat[,1], dev=1:ncol(tri))
    tri <- as.triangle(tri)
    plot(tri)
    plot(tri, lattice=TRUE)

    myCSVfile1 <- "C:/Users/vikto/Desktop/Bakalaurinis/MOD_values.csv" #realių reikšmių po trikampiu nuskaitymas iš failo
    data <- read.csv(file=myCSVfile1)
    triangle <- data[,-1]
    dimnames(triangle) <- list(origin=data[,1], dev=1:ncol(triangle))
    full.triangle <- tri
    n <- ncol(tri)

    f_1 <- sapply((n-1):1, function(i){ #žalų vystymosi koeficiento OLS skaičiavimai
    sum(tri[1:i, n-i+1]*tri[1:i, n-i])/sum(tri[1:i, n-i]^2)
    })
    f_1 <- replace(f_1, f_1<=1, 1.0001)
    dev_period <- 1:(n-1)
    plot(log(f_1-1) ~ dev_period, main="Log-linear extrapolation of age-to-age factors")
    tail.model <- lm(log(f_1-1) ~ dev_period)
    abline(tail.model)
    co <- coef(tail.model)
    tail <- exp(co[1] + c(n:(n + 100)) * co[2]) + 1 #ekstrapoliavimas 100 vystymosi periodų
    f_1.tail <- prod(tail)
    f_1.tail
    (f_1 <- c(f_1, f_1.tail))

    plot(100*(rev(1/cumprod(rev(c(f_1, tail[tail>1.0001]))))), t="b",
    main="Expected claims development pattern",
    xlab="Dev. period", ylab="Development % of ultimate loss") #grafikas

    f_2 <- sapply((n-1):1, function(i){ #žalų vystymosi koeficiento CL skaičiavimai
    sum(tri[1:i, n-i+1])/sum(tri[1:i, n-i])
    })
    f_2 <- replace(f_2, f_2<=1, 1.0001)
    dev_period <- 1:(n-1)
    plot(log(f_2-1) ~ dev_period, main="Log-linear extrapolation of age-to-age factors")
    tail.model <- lm(log(f_2-1) ~ dev_period)
    abline(tail.model)
    co <- coef(tail.model)
    tail <- exp(co[1] + c(n:(n + 100)) * co[2]) + 1 #ekstrapoliavimas 100 vystymosi periodų
    f_2.tail <- prod(tail)
    f_2.tail
    (f_2 <- c(f_2, f_2.tail))

    f <- sapply((n-1):1, function(i){ #žalų vystymosi koeficiento ADF skaičiavimai
    tri[1:i, n-i+1]/tri[1:i, n-i]
    })
    f <- tapply(unlist(f), names(unlist(f)), sum)
    f <- f[2:n]
    f <- unname(f)
    f_3 <- c()
    for(i in 1:n) {
    d <- n-i+1
    f_3 <- c(f_3, (1/d)*f[i])
    }
    f_3 <- f_3[1:(n-1)]
    f_3 <- replace(f_3, f_3<=1, 1.0001)
    (f_3 <- c(f_3, 1))


    simulated.triangle <- full.triangle #sukuriamas naujas trikampis rezultatams talpinti
    simulated.triangle[upper.tri(simulated.triangle)] <- 0
    simulated.triangle[lower.tri(simulated.triangle)] <- 0
    simulated.triangle[row(simulated.triangle) == col(simulated.triangle)] <-0

    func <- function(pe, h) { #skaičiavimai pagal modelį su pasirinktais parametrais
    n <- ncol(tri)
    for(k in 1:(n-1)){
    full.triangle[(n-k+1):n, k+1] <- round(full.triangle[(n-k+1):n, k]*f_1[k]+full.triangle[(n-k+1):n,k]*h*pe, digits=2)
    simulated.triangle[(n-k+1):n, k+1]<- full.triangle[(n-k+1):n, k+1] 
    }

    #modelis A
    #func <- function(pe) { #skaičiavimai pagal modelį A su pasirinktais parametrais 
    #n <- ncol(tri)
    #alpha <- 3
    #sigma <- 1
    #for(k in 1:(n-1)){
    #full.triangle[(n-k+1):n, k+1] <- round(full.triangle[(n-k+1):n, k]*f_1[k]+sigma*(full.triangle[(n-k+1):n,k]^(1-alpha/2))*pe, digits=2)
    #simulated.triangle[(n-k+1):n, k+1]<- full.triangle[(n-k+1):n, k+1] 
    #}

    #modelis B
    #func <- function(pe) { #skaičiavimai pagal modelį B su pasirinktais parametrais 
    #n <- ncol(tri)
    #a <- 3
    #sigma <- 0.01
    #for(k in 1:(n-1)){
    #full.triangle[(n-k+1):n, k+1] <- round(full.triangle[(n-k+1):n, k]*f_1[k]+sigma*sqrt(full.triangle[(n-k+1):n,k]+a*full.triangle[(n-k+1):n,k]^2)*pe, digits=2)
    #simulated.triangle[(n-k+1):n, k+1]<- full.triangle[(n-k+1):n, k+1] 
    #}


    l <- list() #sukuriamas naujas sąrašas
    for(row in 1:nrow(simulated.triangle)) { #reikšmės sudedamos į sąrašą
        for(col in 1:ncol(simulated.triangle)) {
            l <- append(l, simulated.triangle[row, col])
    }}
    return(list("v_1"=l[1], #grąžinamos reikšmės
    "v_2"=l[2],
    ...
    "v_1444"=l[1444]
    ))
    }

    pe <- runif(500, min = -sqrt(3), max = sqrt(3)) #a.d. epsilon generavimas iš tolygiojo skirstinio
    pe <- (pe-mean(pe))/sqrt(var(pe)) 

    a <- 0 #parametrai pagrindiniam modeliui
    alpha <- 0
    sigma <- 0.01
    x <- 1
    h <- sigma*x*sqrt(x^(-alpha)+a)
    param_list <- list("pe"=pe, "h"=h) #parametrų sąrašas pagrindiniam modeliui
    #param_list <- list("pe"=pe) #parametrų sąrašas modeliams A arba B
    erg <- MonteCarlo(func=func, nrep=1, param_list=param_list, max_grid = 10000) #Monte Karlo funkcija

    head(MakeFrame(erg)) 
    df <- MakeFrame(erg)
    df <- round(colMeans(df),digits=2)
    df <- data.matrix(df, rownames.force = NA)
    to_columns <- df[3:nrow(df)]
    #to_columns <- df[2:nrow(df)] naudojamas modeliams A arba B, nes kintamųjų skaičius mažesnis - nėra h funkcijos

    mean_triangle <- matrix(to_columns, nrow=n, ncol = n, byrow = TRUE) #suvidurkintas rezultatų trikampis

    substraction <- round((triangle - mean_triangle), digits=2)^2 #istorinių duomenų ir simuliuotų duomenų vidurkių trikampio atimtis, pakelta kvadratu ir suapvalinta 2 skaičiais po kablelio
    substraction[is.na(substraction)] <-0 #NA reikšmėms priskiriami nuliai
    sums <- colSums(substraction != 0) #jei atimtis nelygi nuliui, tai sumuojama, kiek tokių narių yra
    n <- sum(sums)

    s <- sqrt(1/n * sum(substraction)) #standartinis nuokrypis
    miu <- sum(mean_triangle)/n
    CV <- s/miu*100 #variacijos koeficientas
