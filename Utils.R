# rimuovi caratteri speciali che possono essere problematici in caso di join
clean.text <- function(df,colname,er){
    colname <- enquo(colname)
    df %>% mutate( !!colname := map_chr(name, ~ gsub(er, "", . ))) %>%
        mutate( !!colname := tolower(!!colname)) %>%
        mutate( !!colname := map_chr(name, ~ gsub("\\s+", " ", . )))
}


# restituisce il numero di classi basate sul raggruppamento indicato dalle colonne in ...
number.of.classes <- function(df, ...) nrow(df %>% group_by(...) %>% count())

# conta il numero di linee valide
count.selected.lines <- function(df, ...){
    nrow(df %>% filter(...))
}

# trasforma un campo (trasformabile in numerico) in una serie di fattori (suddivisioni indicate da breaks)
factorise <- function( df, colname, breaks ){
    colname <- enquo(colname)
    df %>% mutate(!!colname := as.numeric(!!colname)) %>%
        mutate( !!colname := cut(!!colname, breaks=breaks) )
}

# filtra per match di espressione regolare
string.query <- function(df, colname, er){
    colname <- enquo(colname)
    df %>% filter(grepl(er, !!colname))
}

# prendi il df degli elementi unici di una colonna
get.unique <- function(df, colname){
    colname <- enquo(colname)
    df %>% select(!!colname) %>% unique()
}

# filter applicato su colonne con tag
filter.by.tag.or <- function(df, colname, ...){
    colname <- enquo(colname)
    # usa degli id per selezionare gli out validi
    df2 <- df %>% mutate(id.. = seq(1,nrow(df))) %>% 
        unnest(!!colname) %>% filter(!!colname %in% ...) %>%
        select(id..) %>% unique()
    semi_join( df %>% mutate(id.. = seq(1,nrow(df))),
               df2, by=c("id..")) %>%
        select(-id..)
}

# filter applicato su colonne con tag
filter.by.tag.negated.or <- function(df, colname, ...){
    colname <- enquo(colname)
    # usa degli id per selezionare gli out validi
    df2 <- df %>% mutate(id.. = seq(1,nrow(df))) %>% 
        unnest(!!colname) %>% filter((!!colname %in% ...)) %>%
        select(id..) %>% unique()
    anti_join( df %>% mutate(id.. = seq(1,nrow(df))),
               df2, by=c("id..")) %>%
        select(-id..)
}

# come sopra ma devono esserci tutte le tag
filter.by.tag.and <- function(df, colname, ...){
    colname <- enquo(colname)
    for(x in ...) filter(df, x %in% !!colname)
}

# funzione per verificare se un grafo è regolarizzabile 
# (slide lezione)
regularify = function (g) {
    n = vcount(g)
    m = ecount(g)
    E = get.edges(g, E(g))
    B = matrix(0, nrow = n, ncol = m)
    # build incidence matrix
    for (i in 1:m) {
        B[E[i,1], i] = 1
        B[E[i,2], i] = 1
    }  
    # objective function
    obj = rep(0, m + 1)
    # constraint matrix
    con = cbind(B, rep(-1, n))
    # direction of constraints
    dir = rep("=", n)
    # right hand side terms
    rhs = -degree(g)
    # solve the LP problem
    sol = lp("max", obj, con, dir, rhs)
    # get solution
    if (sol$status == 0) {
        s = sol$solution
        # weights
        w = s[1:m] + 1
        # weighted degree
        d = s[m+1]
    }
    # return the solution
    if (sol$status == 0) return(list(weights = w, degree = d)) else return(NULL)   
}


# Compute power x = (1/x) A 
#INPUT
# A = graph adjacency matrix
# t = precision
# OUTPUT
# A list with:
# vector = power vector
# iter = number of iterations
# (slide lezione)
power_utils = function(A, t) {
    n = dim(A)[1];
    # x_2k
    x0 = rep(0, n);
    # x_2k+1
    x1 = rep(1, n);
    # x_2k+2
    x2 = rep(1, n);
    diff = 1
    eps = 1/10^t;
    iter = 0;
    while (diff > eps) {
        x0 = x1;
        x1 = x2;
        x2 = (1/x2) %*% A;
        diff = sum(abs(x2 - x0));
        iter = iter + 1;
    } 
    # it holds now: alpha x2 = (1/x2) A
    alpha = ((1/x2) %*% A[,1]) / x2[1];
    # hence sqrt(alpha) * x2 = (1/(sqrt(alpha) * x2)) A
    x2 = sqrt(alpha) %*% x2;
    return(list(vector = as.vector(x2), iter = iter))
}  

# funzione per il calcolo della similarità
# (slide lezione)
similarity = function(g, type = "cosine", mode = "col" ) {
    A = as_adjacency_matrix(g, sparse = FALSE)
    if (mode == "row") {A = t(A)}
    if (type == "cosine") {
        euclidean = function(x) {sqrt(x %*% x)}
        d = apply(A, 2, euclidean)
        D = diag(1/d)
        S = D %*% t(A) %*% A %*% D
    }
    if (type == "pearson") {
        S = cor(A)
    }
    return(S)
}

# Esegue un passo di una visita BFS
BFS.reachable <- function(current, A){
    for(i in current){
        if(i!=0){
            current <- as.numeric(current | A[i,])
        }
    }
    current
}

# versione iterativa della procedura per il calcolo della similarità per tag
exp.corr.similarity.iter <- function(A, selected, SIM_M, max_iter, current, COV){
    done <- rep(0, times=length(current))
    SIM_M[selected,selected] <- 1
    for(iter in 1:max_iter){
        if(iter != 1)
            done <- current
        current <- BFS.reachable(current,A)
        for(j in 1:length(current)){
            if( (current[j] >= 1) && (done[j] == 0)){
                SIM_M[selected,j] <- max(SIM_M[selected,j], COV[selected,j]/(2^(iter-1)))
                SIM_M[j,selected] <- SIM_M[selected,j]
            }
        }
    }
    return(SIM_M)
}

# versione ricorsiva (più lenta) della procedura per il calcolo della similarità per tag
exp.corr.similarity <- function(A, bit, selected, SIM_M, iter, max_iter, current, COV){
    if(iter == 0){
        SIM_M[selected,selected] <- 1
        return(exp.corr.similarity(A, bit, selected, SIM_M, iter+1, max_iter, current, COV))
    }
    if(iter <= max_iter){
        current <- BFS.reachable(current,A)
        for(j in 1:length(current)){
            if( (current[j] >= 1)){
                SIM_M[selected,j] <- max(SIM_M[selected,j], COV[selected,j]/(2^(iter-1)))
                SIM_M[j,selected] <- SIM_M[selected,j]
            }
        }
        return(exp.corr.similarity(A, bit, selected, SIM_M, iter+1, max_iter, current, COV))
    }else{
        return(SIM_M)
    }
}

# calcolo dell'agreement tra bitset
compute.bitcor <- function(bitv){
    m <- matrix( rep(0,times=nrow(bitv)^2), byrow = TRUE, ncol=nrow(bitv))
    for(i in 1:nrow(m)){
        for(j in i:ncol(m)){
            # calcola l'agreement %
            m[i,j] <- (ncol(bitv) - sum(xor(bitv[i,], bitv[j,])))/ncol(bitv)
            m[j,i] <- m[i,j]
        }
    }
    m
}

# media la somma dei valori di un vettore dividendo 
# per i soli numeri positivi
mean.of.positives <- function(M){
    total <- apply(M, 1, sum)
    positives <- apply(M,1,
                       function(l){
                           cnt <- 0
                           for(i in l)
                               if(i>0) cnt <- cnt + 1
                               ifelse(cnt>0,cnt,1)
                       })
    total/positives
}


# Calcola l'affinità media di un nodo rispetto alle conponenti di una selezione
# sfruttando la compatibilità di bitsets associati ai nodi
#
# g : grafo G(V,E)
# sel : nodi di V da scansionare
# COV : matrice precalcolata delle compatibilità dei bitset
# depth : profondità della BFS per la propagazione dell'affinità con un singolo nodo  
#
get.affinity <- function(g, sel, COV, depth){
    A <- as_adjacency_matrix(g, sparse = FALSE)
    Sim_matrix <- matrix(rep(0, times=nrow(A)^2), ncol=nrow(A))
    for(i in sel){
        ini <- rep(0,times=ncol(A))
        ini[i] <- 1
        Sim_matrix <- exp.corr.similarity.iter(A, i, Sim_matrix, depth, ini, COV)
    }
    mean.of.positives(Sim_matrix)
}

#### shared Userbase

# calcola la lista di adiacenza dalla matrice
to.adj.list<-function(A){
    AL <- list()
    for(i in 1:nrow(A)){
        AL[[i]]<-which(A[i,] == 1)
    }
    AL    
}

# calcolo della similarità per userbase condivisa (versione R)
get.shared.userbase.similarity <- function(g, sel, depth, attr){
    A <- as_adjacency_matrix(g, sparse = FALSE)
    AL <- to.adj.list(A)
    Sim_matrix <- matrix(rep(0, times=nrow(A)*length(sel)), ncol=nrow(A))
    x <- 1
    for(i in sel){
        ini <- set(i)
        Sim_matrix[x] <- shared.userbase.similarity(A,AL, i, depth, ini, g, attr)
        x <- x+1
    }
    apply(Sim_matrix, 2, sum)
}

# funzione ausiliaria per il calcolo della similarità sulla base di uno specifico gioco/nodo
shared.userbase.similarity <- function(A,AL, i, depth, Q, g, attr){
    flow <- rep(0,times=ncol(A))
    done <- rep(FALSE,times=ncol(A))
    while(!set_is_empty(Q)){
        # prendi l'elemento dalla coda
        for(v in Q){
            print(as.list(v)[[1]])
            Q <- Q-v
            break
        }
        
        for(d in AL[[v]]){
            if(!done[d]){
                w <- get.edge.attribute(g,attr,
                       get.edge.ids(g,c(v,d)))
                flow[d] <- w + flow[v]/2
                Q<-Q+d
            }
        }
        done[v]<-TRUE
    }
    return(flow)
}

# funzione per il calcolo della cosine similarity media
# per un insieme di nodi e un altro nodo
cosine.set.similarity <- function(sel, cosine.similarity.matrix){
    cosine <- rep(0,times=ncol(cosine.similarity.matrix))
    for(i in 1:length(cosine)){
        cosine[i] <- sum(map_dbl(sel, ~cosine.similarity.matrix[i,.]))/length(sel)
    }
    cosine
}

