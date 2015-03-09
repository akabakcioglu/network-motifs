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



# --- calculate average cooperativity over random networks
n = 11
ntotal = 100000;            # ensemble size
numinit = 2^11              # num of initial cond for fp determination

for ne=29:29
    for np=15:15
        dirname = "$n.$ne.$np"
        if (!isdir(dirname))
            run(`mkdir $dirname`)
        end
        nself = 5
        selfsh = [ones(Int64,nself),zeros(Int64,n-nself)]
        listc=zeros(Float64,ntotal)        # list of coherence values
        lista=zeros(Float64,ntotal)        # list of active node fractions
        listl=zeros(Float64,ntotal)        # list of fraction of point attractors
        intersh = randomnkp(selfsh,ne,np) # generate one sample with given (nkp)
        i=0
        while (i<ntotal)
            (intersh,selfsh) = shufflenet(intersh,selfsh) # shuffle
            if (isconnected(symadj(intersh))) # check if connected
#                (ncop,actratio,natt,npointattr) = coopattr(intersh,selfsh,numinit)
                (attlist,basinlist,listperiod,listcop,listact) = coopattr(intersh,selfsh,numinit)
                if (length(attlist) > 1)
                    i += 1              # increase network count
                    edges = zeros(Int64,ne+sum(selfsh),3)
                    idx=1
                    for j=1:n
                        if selfsh[j] == 1
                            edges[idx,:] = [j j -1] # self degrading genes
                            idx += 1
                        end
                        for k=1:n
                            if intersh[j,k] != 0
                                edges[idx,:] = [j k intersh[j,k]]
                                idx += 1
                            end end end
                    writedlm("$dirname/$i.dat",edges,' ')
                    writedlm("$dirname/$i.att",[attlist' basinlist listperiod listcop/n listact],' ')
                    if (mean(listcop/n)==1)
                        print("$i\n")
                    end
                end
            end
        end
    end
end
