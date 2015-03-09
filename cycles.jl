srand(1)

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



function step (nw,pth,visnodes,cdict)

    sizenw = length(nw[:,1])
    p = copy(pth)
    last = pop!(p)
    if in(last,p)                       # found cycle
        pp = copy(p)                    # make a copy of the path (something like 1-2-10-3)
        c = Int64[]       # cycle nodes ordered with correct direction
        mini = last       # stores the minimum index in cycle
        unshift!(c,last)
        n=pop!(pp)
        while n != last
            unshift!(c,n)
            if(n<mini)
                mini = n
            end
            n=pop!(pp)
        end

        while c[1] != mini      # bring the min index to the beginning
            n = shift!(c)
            push!(c,n)
        end
        cyckey = join(c,"-")

        lcyc = length(c)
        if (!haskey(cdict,cyckey))
            cdict[cyckey]=1
            for i=1:lcyc-1
                cdict[cyckey] *= nw[c[i+1],c[i]]
            end
            cdict[cyckey] *= nw[c[1],c[lcyc]]
        end
    else                                # move one more step if possible
        push!(p,last)
        for i=1:sizenw
            if nw[i,last] != 0
                push!(p,i)
                visnodes[i] = 1
                step(nw,p,visnodes,cdict)
                pop!(p)
            end
        end
    end
end


function getcycles(nw)
    visnodes = Dict()                       # set of visited nodes
    path=Int64[]                            # ordered list of path nodes
    cycdict = Dict()                        # set of cycles found so far
    sizenw=length(nw[:,1])

    ncycles=0
    for j=1:sizenw
        if !haskey(visnodes,j)
            push!(path,j)
            step(nw,path,visnodes,cycdict)
            pop!(path)
        end
    end
    return cycdict
end


nw = readdlm("7org_nw/d.melanogaster")      # interaction matrix
cycdict = getcycles(nw)
println(sum(values(cycdict))/length(cycdict),"\t",length(cycdict))

## # --- calculate average cooperativity over random networks
n = 0
ntotal = 1000000           # ensemble size
numinit = 2^11              # num of initial cond for fp determination

path="7org_nw"
## ## ## for f=["a.thaliana","d.melanogaster","lymphocyte","mammalian","myeloid","s.cerevisiae","s.pombe"
#for f=["a.thaliana","lymphocyte","mammalian","myeloid","s.cerevisiae","s.pombe"]
for f=["s.cerevisiae","s.pombe"]
    nw = readdlm("$path/$f")
    println("&$f")
    selfsh = diag(nw)
    intersh = nw-diagm(selfsh)
    n = length(nw[:,1])
    i=0
    while (i<ntotal)
        (intersh,selfsh) = shufflenet(intersh,selfsh) # shuffle
        if (isconnected(symadj(intersh))) # check if connected
            (attlist,basinlist,listperiod,listcop,listact) = coopattr(intersh,selfsh,numinit)
            if (length(attlist) > 1)
                i += 1              # increase network count
                nw = intersh + diagm(selfsh)
                cycdict = getcycles(nw)
                cyccount = length(cycdict)
                cycscore = sum(values(cycdict))/cyccount
                println("$i\t$cyccount\t$cycscore")
                if (mean(listcop/n)==1) # perfect coherence
                    println("*$i\t$cyccount\t$cycscore")
                end
            end
        end
    end
end
