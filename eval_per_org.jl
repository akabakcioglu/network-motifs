srand(1)
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

# generate motif id maps
(motifid,motifidrev) = motifmap();
(motifidnoloop,motifidnolooprev) = motifmapnoloop();

path="7org_nw"
f=pop!(ARGS)
nw = int(readdlm("$path/$f"))
self = diag(nw)
selfm = diagm(self)
inter = nw-selfm
n = length(inter[:,1])
ne = sum(inter .* inter)
np = int((ne+sum(inter))/2)

dirname = "$f.$n.$ne.$np"
motifcount = readdlm("$dirname/loop.$f")
motifcountnoloop = readdlm("$dirname/noloop.$f")
motifcountperf = readdlm("$dirname/perf.$f")
motifcountnoloopperf = readdlm("$dirname/perfnl.$f")
motifcountbio = readdlm("$dirname/bio.$f")
motifcountnoloopbio = readdlm("$dirname/bionl.$f")

motifcount /= sum(motifcount)
motifcountnoloop /= sum(motifcountnoloop)
motifcountperf /= sum(motifcountperf)
motifcountnoloopperf /= sum(motifcountnoloopperf)
motifcountbio /= sum(motifcountbio)
motifcountnoloopbio /= sum(motifcountnoloopbio)
logp_bio_all=0
for i=1:length(motifcount)
    if (motifcountbio[i]>0 && motifcount[i]>0)
        logp_bio_all += motifcountbio[i]*log(motifcount[i])
    else
        if (motifcountbio[i]>0 && motifcount[i]==0)
            print("bad motif in bio: $i: ")
            println(base(3,motifidrev[i],9))
        end
    end
end
println("log likelihood for $f in null set = $logp_bio_all")
logp_bio_perf=0
for i=1:length(motifcountperf)
    if (motifcountbio[i]>0 && motifcountperf[i]>0)
        logp_bio_perf += motifcountbio[i]*log(motifcountperf[i])
    else
        if (motifcountbio[i]>0 && motifcountperf[i]==0)
            print("bad motif in bio (perf): $i: ")
            println(base(3,motifidrev[i],9))
        end
    end
end
println("log likelihood for $f in perfect set = $logp_bio_perf")
logp_bio_all_nl=0
for i=1:length(motifcountnoloop)
    if (motifcountnoloopbio[i]>0 && motifcountnoloop[i]>0)
        logp_bio_all_nl += motifcountnoloopbio[i]*log(motifcountnoloop[i])
    else
        if (motifcountnoloopbio[i]>0 && motifcountnoloop[i]==0)
            print("bad motif in bio noloop: $i: ")
            println(base(3,motifidrev[i],9))
        end
    end
end
println("log likelihood for $f in null set (noloop) = $logp_bio_all_nl")
logp_bio_perf_nl=0
for i=1:length(motifcountnoloopperf)
    if (motifcountnoloopbio[i]>0 && motifcountnoloopperf[i]>0)
        logp_bio_perf_nl += motifcountnoloopbio[i]*log(motifcountnoloopperf[i])
    else
        if (motifcountnoloopbio[i]>0 && motifcountnoloopperf[i]==0)
            print("bad motif in bio (noloop-perf): $i: ")
            println(base(3,motifidrev[i],9))
        end
    end
end
println("log likelihood for $f in perfect set no loop = $logp_bio_perf_nl")
logp_perf_all=0
for i=1:length(motifcount)
    if (motifcountperf[i]>0 && motifcount[i]>0)
        logp_perf_all += motifcountperf[i]*log(motifcount[i])
    else
        if (motifcountperf[i]>0 && motifcount[i]==0)
            print("bad motif in perf: $i: ")
            println(base(3,motifidrev[i],9))
        end
    end
end
println("log likelihood for $f-perf in null set = $logp_perf_all")
logp_perf_all_nl=0
for i=1:length(motifcountnoloop)
    if (motifcountnoloopperf[i]>0 && motifcountnoloop[i]>0)
        logp_perf_all_nl += motifcountnoloopperf[i]*log(motifcountnoloop[i])
    else
        if (motifcountnoloopperf[i]>0 && motifcountnoloop[i]==0)
            print("bad motif in perf noloop: $i: ")
            println(base(3,motifidrev[i],9))
        end
    end
end
println("log likelihood for $f-perf in null set (noloop) = $logp_perf_all_nl")
