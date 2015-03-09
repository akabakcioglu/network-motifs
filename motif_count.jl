srand(1)


# --- generate a random (n,ne,np) network.
# --- d = matrix diagonal
# --- ne = number of edges
# --- np = number of positive edges
function randomnkp (d,ne,np)
    n = length(d)                       # int matrix diagonal - never change
    nn = ne-np                          # negative edges
    nz = n*(n-1)-ne                     # number of zero-edges
    intlist = [ones(Int64,np);-ones(Int64,nn);zeros(Int64,nz)]
    rect = reshape(shuffle(intlist),n,n-1)
    rrect = [rect zeros(Int64,n)]       # now insert back the diagonal
    lrect = [zeros(Int64,n) rect]
    interrand = tril(rrect - diagm(diag(rrect))) + triu(lrect-diagm(diag(lrect)))
    return interrand
end


# --- return a (n,ne,np) preserved shuffled version of a network
function shufflenet (intr,slf)
    n = length(slf)
    nw = similar(intr)                  # this will be the randomized version
    shf = shuffle(slf)                  # keep number of self degredating genes
    tl = tril(intr)                     # remove diag and join upper/lower triang
    tu = triu(intr)
    rect = tu[:,2:end]+tl[:,1:end-1]    # shuffle this
    rect = reshape(shuffle(reshape(rect,length(rect),1)[:,1]),n,n-1)
    rrect = [rect zeros(Int64,n)]       # now insert back the shuffled diagonal
    lrect = [zeros(Int64,n) rect]
    nw = tril(rrect - diagm(diag(rrect))) + triu(lrect-diagm(diag(lrect)))
    return (nw,shf)
end


# --- generate the symmetric adjacency matrix with elements 0,1 ---
function symadj(m)
    n = length(m[:,1])
    ad = zeros(Int64,n,n)
    for (i=1:n) for (j=i+1:n)       # adjacency matrix - without diagonals
        ad[i,j] = m[i,j]==0 && m[j,i]==0 ? 0 : 1
        ad[j,i] = ad[i,j]
    end end
    return ad
end


## --- check connectivity -----------
function isconnected (adj)
    n = length(adj[:,1])                # number of nodes
    nsq = n*n
    dum = adj
    (sumprev,sumdum) = (0,0)
    for (i=1:n)
        dum = sign(dum + dum*adj)
        sumdum = sum(dum)
        if (sumdum == sumprev)
            break
        end
        sumprev = sumdum
    end
    return (sumdum==nsq)
end    


## --- node update function (integrator): yeast's majority rule -----------
function upd(s,intr,slf,n)
    sss = similar(s)
    inp = intr*s
    for (i=1:n)
        sss[i] = (inp[i]>0)?1:(inp[i]<0)?0:(slf[i]==1)?0:s[i]
    end
    return sss
end


## ---- output is the attractor list + basin sizes ---
## ---- calculated from a subset of initial conditions ---
function getcyclessample(intr,slf,ninit)
    
    n = length(slf)
    attrlist = zeros(Int64,n)           # s=0 is always an attractor
    keyzeros = join(zeros(Int64,n))
    attdict = [keyzeros=>0]    # values are the basin sizes
    for (i=1:ninit)       # calculate basins using ninit initial states
        s = rand(0:1,n)
        skey = join(s)
        pathdict = Dict()
        while (!haskey(attdict,skey) && !haskey(pathdict,skey))
            pathdict[skey] = "-"
            s = upd(s,intr,slf,n)
            skey=join(s)
        end
        if (haskey(attdict,skey)) attdict[skey] += 1 
        else                        # this is a new cycle or fp
            attdict[skey] = 1
            attrlist = [attrlist s]
        end
    end
    basinlist = [attdict[join(attrlist[:,i])]/ninit for i in 1:length(attrlist[1,:])]
    return attrlist,basinlist
end


## ---- find number of coop nodes in a state ---
function ncoop (s,intr)
    n = length(s)
    act = intr*diagm(s) # keep only the columns for active nodes
    nc=0
    for (i=1:n)
#        if (abs(sum(act[i,:])) == sum(abs(act[i,:])))
        if (!in(1,act[i,:]) || !in(-1,act[i,:]))
            nc += 1
        end
    end
    return nc
end


## ---- find coop coef for sampled cycles of a network ---
function coopattr (intr,slf,ninit)
    n = length(slf)
    (listatt,listb) = getcyclessample(intr,slf,ninit)
    listperiod=Int64[];
    listcop=Float64[];
    listact=Float64[];
    npointattr = 0;
    for (i=1:length(listatt[1,:]))        # go over each attractor
        s = listatt[:,i]
#        if (sum(s) != 0)
        ncp = ncoop(s,intr)
#            nac = sum(abs(intr)*s)
        nac = sum(s)
        period = 1
        s=upd(s,intr,slf,n)
        while (s != listatt[:,i])
            ncp += ncoop(s,intr)
#            nac += sum(abs(intr)*s)
            nac += sum(s)
            period += 1
            s=upd(s,intr,slf,n)
        end
        push!(listperiod,period)
        push!(listcop,ncp/period)
        push!(listact,nac/period)
        if (period==1) npointattr += 1 end
    #        end
    end
#    if (length(listperiod)==1)  return -1,-1,-1,-1 # no fixed point other than (0,..0)
#    else                    # also returns the fracion of active nodes
#        return (1.0*mean(listcop)/n,mean(listact)/n,length(listc[1,:]),npointattr)
    return listatt,listb,listperiod,listcop,listact
#    end
end

# returns the node-id permuted version of the square matrix a
function motifpermute (a,perm)
    b = similar(a)
    n = length(a[:,1])
    for i=1:n
        for j=1:n
            b[i,j]=a[perm[i],perm[j]]
        end
    end
    return b
end


# returns the dumb id of a motif matrix (3x3)
function dumbid (a)
    id=0;
    for i=1:3
        for j=1:3
            digit = 3*(i-1)+(j-1) # (1,1)=0, (1,2)=1,..,(2,1)=3..
            id += (a[i,j]+1)*(3^digit)
        end
    end
    return id
end


function motifmap ()
    a=zeros(Int64,3,3)
    motifid = Dict()                    # values are the true ids
    motifidrev = Dict()                 # values are the dumb ids
    idx = 1 # this will be used as an array index, so I don't start from 0.
    for a[1,1]=-1:1
        for a[1,2]=-1:1
            for a[1,3]=-1:1
                for a[2,1]=-1:1
                    for a[2,2]=-1:1
                        for a[2,3]=-1:1
                            for a[3,1]=-1:1
                                for a[3,2]=-1:1
                                    for a[3,3]=-1:1
                                        perm = [1,2,3]
                                        id = dumbid(motifpermute(a,perm))
                                        if (!haskey(motifid,id))
                                            motifidrev[idx]=id
                                            for perm in collect(permutations([1:3]))
                                                motifid[dumbid(motifpermute(a,perm))]=idx
                                            end
                                            idx += 1
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return motifid,motifidrev
end


                                    
function motifmapnoloop ()
    a=zeros(Int64,3,3)
    motifid = Dict()                    # values are the true ids
    motifidrev = Dict()                 # values are the dumb ids
    idx=1 # this will be used as an array index, so I don't start from 0.
    for a[1,2]=-1:1
        for a[1,3]=-1:1
            for a[2,1]=-1:1
                for a[2,3]=-1:1
                    for a[3,1]=-1:1
                        for a[3,2]=-1:1
                            perm = [1,2,3]
                            id = dumbid(motifpermute(a,perm))
                            if (!haskey(motifid,id))
                                motifidrev[idx]=id
                                for perm in collect(permutations([1:3]))
                                    motifid[dumbid(motifpermute(a,perm))]=idx
                                end
                                idx += 1
                            end
                        end
                    end
                end
            end
        end
    end
    return motifid,motifidrev
end
                                    

function getmotifstats(nw,map,mc)
    
    n=length(nw[:,1])
    a=zeros(Int64,3,3)
    for i=1:n
        for j=i+1:n
            for k=j+1:n
                a = nw[[i,j,k],[i,j,k]]
                mc[map[dumbid(a)]] += 1
            end
        end
    end
end

# --- calculate average cooperativity over random networks
n = 11
ntotal = 1000000;           # ensemble size
numinit = 2^11              # num of initial cond for fp determination

# generate motif id maps
(motifid,motifidrev) = motifmap();

a = zeros(Int64,0)
for x in values(motifid) push!(a,x) end
## motifcount=zeros(Int64,sort(a)[length(a)])
## motifcountperf=zeros(Int64,sort(a)[length(a)])
## motifcountbio=zeros(Int64,sort(a)[length(a)])

(motifidnoloop,motifidnolooprev) = motifmapnoloop();

a = zeros(Int64,0)
for x in values(motifidnoloop) push!(a,x) end
## motifcountnoloop=zeros(Int64,sort(a)[length(a)])
## motifcountnoloopperf=zeros(Int64,sort(a)[length(a)])
## motifcountnoloopbio=zeros(Int64,sort(a)[length(a)])

## ## for f=["a.thaliana","d.melanogaster","lymphocyte","mammalian","myeloid","s.cerevisiae","s.pombe"]
## for f=["a.thaliana","lymphocyte","mammalian","myeloid","s.cerevisiae","s.pombe"]
##     path="7org_nw"
##     nw = readdlm("$path/$f");
##     getmotifstats(nw,motifid,motifcountbio)
##     intersh = nw-diagm(diag(nw));
##     getmotifstats(intersh,motifidnoloop,motifcountnoloopbio)
## end
## for ne=29:29
##     for np=15:15
##         ## dirname = "$n.$ne.$np"
##         ## if (!isdir(dirname))
##         ##     run(`mkdir $dirname`)
##         ## end
##         nself = 5
##         selfsh = [ones(Int64,nself),zeros(Int64,n-nself)]
##         listc=zeros(Float64,ntotal)    # list of coherence values
##         lista=zeros(Float64,ntotal)    # list of active node fractions
##         listl=zeros(Float64,ntotal) # list of fraction of point attractors
##         intersh = randomnkp(selfsh,ne,np) # generate one sample with given (nkp)
##         i=0
##         while (i<ntotal)
##             (intersh,selfsh) = shufflenet(intersh,selfsh) # shuffle
##             if (isconnected(symadj(intersh))) # check if connected
##                 (attlist,basinlist,listperiod,listcop,listact) = coopattr(intersh,selfsh,numinit)
##                 if (length(attlist) > 1)
##                     i += 1              # increase network count
##                     nw = intersh - diagm(selfsh)
##                     ## edges = zeros(Int64,ne+sum(selfsh),3)
##                     ## idx=1
##                     ## for j=1:n
##                     ##     for k=1:n
##                     ##         if nw[j,k] != 0
##                     ##             edges[idx,:] = [j k nw[j,k]]
##                     ##             idx += 1
##                     ##         end end end
##                     ## writedlm("$dirname/$i.dat",edges,' ')
##                     ## writedlm("$dirname/$i.att",[attlist' basinlist listperiod listcop/n listact],' ')
##                     getmotifstats(nw,motifid,motifcount)
##                     getmotifstats(intersh,motifidnoloop,motifcountnoloop)
##                     if (mean(listcop/n)==1) # perfect coherence
##                         print("$i\n")
##                         getmotifstats(nw,motifid,motifcountperf)
##                         getmotifstats(intersh,motifidnoloop,motifcountnoloopperf)
##                     end
##                 end
##             end
##         end
##     end
## end
