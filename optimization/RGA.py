
from optimization.calcFunc import evaluateAero 
import numpy as np
import matplotlib.pyplot as plt

def soga(fun_name,nvar,ncon,lb,ub,npop,maxgen,pcross,pmut,cht_type=None):
    # Obtaining initial parent population (normalized form with value between 0 and 1)
    pop_parent = np.random.rand(npop,nvar)

    # Calculating responses (objective and constraint values) of parent population
    objval_parent,conval_parent = calc_obj_con(pop_parent,fun_name,ncon,lb,ub,cht_type)

    # Defining fitness value (can be real values, rank from MCR, etc) of parent poopulation
    fitval_parent,cv_parent,Nv_parent = evaluate_fitval(objval_parent,conval_parent,cht_type)

    # Finding feasible individuals of parent population
    ind_parent_feas = np.where(Nv_parent==0)[0]

    # Allocating lists for current best individual
    best_individual = []
    best_objval = []
    best_conval = []
    best_fitval = []

    if len(ind_parent_feas) > 0:
        pop_parent_feas = pop_parent[ind_parent_feas,:]
        objval_parent_feas = objval_parent[ind_parent_feas]
        conval_parent_feas = conval_parent[ind_parent_feas,:]
        fitval_parent_feas = fitval_parent[ind_parent_feas]

        # Sorting feasible individuals of parent population from lowest to highest fitness values
        allIndex = np.argsort(fitval_parent_feas)

        # Recording the current best individual
        best_individual.append(pop_parent_feas[allIndex[0],:])
        best_objval.append(objval_parent_feas[allIndex[0]])
        best_conval.append(conval_parent_feas[allIndex[0],:])
        best_fitval.append(fitval_parent_feas[allIndex[0]])
    else:
        # Sorting feasible individuals of parent population from lowest to highest fitness values
        allIndex = np.argsort(fitval_parent)

        # Recording the current best individual
        best_individual.append(np.nan*np.ones(nvar))
        best_objval.append(np.nan)
        best_conval.append(np.nan*np.ones(ncon))
        best_fitval.append(np.nan)

    # Generational looping
    for gen in range(maxgen):
        # if gen == 0:
        #     print('currently evaluating gen:',gen,', current optimum: ',objval_parent_feas[allIndex[0]])
        # else:
        #     print('currently evaluating gen:',gen,', current optimum: ',objval_child_feas[allIndex[0]])

        # Stopping criterion
        if gen == maxgen-1:
          break

        # Creating mating pool
        matingpool = create_matingpool(pop_parent,fitval_parent)

        # Conducting crossover
        pop_child = crossover_process(matingpool,pcross)

        # Conducting mutation
        pop_child_mut = mutation_process(pop_child,pmut)

        # Calculating responses (objective and constraint values) of child population
        objval_child,conval_child = calc_obj_con(pop_child_mut,fun_name,ncon,lb,ub,cht_type)

        # Elitism
        if len(ind_parent_feas) > 0:
            pop_child_mut[0,:] = best_individual[len(best_individual)-1]
            objval_child[0] = best_objval[len(best_objval)-1]
            conval_child[0,:] = best_conval[len(best_conval)-1]
        else:
            pop_child_mut[0,:] = pop_parent[allIndex[0],:]
            objval_child[0] = objval_parent[allIndex[0]]
            conval_child[0,:] = conval_parent[allIndex[0],:]

        # Defining fitness value (can be real values, rank from MCR, etc) of combined poopulation
        fitval_child,cv_child,Nv_child = evaluate_fitval(objval_child,conval_child,cht_type)

        # Finding feasible individuals of child population
        ind_child_feas = np.where(Nv_child==0)[0]

        if len(ind_child_feas) > 0:
            pop_child_feas = pop_child[ind_child_feas,:]
            objval_child_feas = objval_child[ind_child_feas]
            conval_child_feas = conval_child[ind_child_feas,:]
            fitval_child_feas = fitval_child[ind_child_feas]

            # Sorting individuals of combined population from lowest to highest fitness values
            allIndex = np.argsort(fitval_child_feas)

            # Recording the current best individual
            best_individual.append(pop_child_feas[allIndex[0],:])
            best_objval.append(objval_child_feas[allIndex[0]])
            best_conval.append(conval_child_feas[allIndex[0],:])
            best_fitval.append(fitval_child_feas[allIndex[0]])
        else:
            # Sorting individuals of combined population from lowest to highest fitness values
            allIndex = np.argsort(fitval_child)

            # Recording the current best individual
            best_individual.append(np.nan*np.ones(nvar))
            best_objval.append(np.nan)
            best_conval.append(np.nan*np.ones(ncon))
            best_fitval.append(np.nan)

        # New population
        pop_parent = pop_child_mut.copy()
        fitval_parent = fitval_child.copy()
        objval_parent = objval_child.copy()
        conval_parent = conval_child.copy()
        Nv_parent = Nv_child.copy()

        # Finding feasible individuals of parent population
        ind_parent_feas = np.where(Nv_parent==0)[0]

        # fig = plt.figure()
        # plt.plot(np.arange(0,gen),best_objval[len(best_objval)-1])
        # plt.xlabel('Generation')
        # plt.ylabel('Obtained Minimum Value')
        # plt.grid()
        # plt.show()

        
        with open(".\objval.txt", "ab") as f:
            np.savetxt(f, best_objval, fmt='  %5f\n')

    best_individual = np.array(best_individual)
    best_objval = np.array(best_objval)
    best_conval = np.array(best_conval)

    result = [best_objval,best_individual,best_conval]
    return result

def calc_obj_con(population,fun_name,ncon,lb,ub,cht_type=None):
    npop,nvar = population.shape
    objval = np.zeros(npop)         # initializing storing array for objective values
    if ncon == 0:
        conval = np.zeros((npop,1)) # dummy constraint for unconstrained problem
    else:
        conval = np.zeros((npop,ncon)) # initializing storing array for constraint values

    for i in range(npop):
        individual = population[i,:]
        indi_denorm = calc_denorm(individual,lb,ub)
        if ncon == 0:
           objval[i] = evaluateAero(indi_denorm,ncon,fun_name)
        else:
           objval[i],conval[i,:] = evaluateAero(indi_denorm,ncon,fun_name)

    return objval,conval

def evaluate_fitval(objval,conval,cht_type=None):
    npop,ncon = conval.shape
    fitval = np.zeros(npop)

    # Calculating constraint violation fron constraint values
    cv = conval.copy()
    cv[cv<0] = 0

    Nv = np.sum(cv>0,axis=1,dtype='int')  # counting number of constraints violated for each individual
    countfeas = np.sum(Nv==0)             # counting the number of feasible individuals

    if cht_type is None:
        fitval[:] = objval
    elif cht_type == 'SoF':
        sumcv = np.sum(cv,axis=1) # sum of constraint violations
        if countfeas > 0:
            objmax = np.amax(objval[np.where(Nv==0)[0]])
        else:
            objmax = 0
        for i in range(npop):
            if sumcv[i] == 0.0:
                fitval[i] = objval[i]
            else:
                fitval[i] = objmax + sumcv[i]
    elif cht_type == 'G-MCR':
        con_check = np.sum(cv>0,axis=0)
        Fobj = np.zeros(npop)
        Fcon = np.zeros((npop,ncon))
        FNv = np.zeros(npop)
        for i in range(npop):
            Fobj[i] = np.sum(objval<objval[i])
            FNv[i] = np.sum(Nv<Nv[i])
            for j in range(ncon):
                Fcon[i,j] = np.sum(cv[:,j]<cv[i,j])

        beta1 = np.sqrt(1-(countfeas/npop-1)**2)
        beta2 = 1-beta1
        eta = 0
        alpha = (1/npop)*con_check
        sum_alpha = np.sum(alpha)
        if sum_alpha > 0:
            gamma = 1/sum_alpha
        else:
            gamma = npop
        fitval[:] = beta1*Fobj + beta2*(eta*FNv + gamma*np.sum(np.multiply(Fcon,alpha),axis=1))

    return fitval,cv,Nv

def create_matingpool(population,fitval):
	npop,nvar = population.shape
	matingpool = np.zeros((npop,nvar))
	for kk in range(npop):
		ip1 = np.random.randint(npop)
		ip2 = np.random.randint(npop)
		if ip2 == ip1:
			while ip2 == ip1:
				ip2 = np.random.randint(npop)
		Ft1 = population[ip1-1,:]
		Ft2 = population[ip2-1,:]
		Fit1 = fitval[ip1-1]
		Fit2 = fitval[ip2-1]
		if Fit1 < Fit2:
			matingpool[kk,:] = Ft1
		else:
			matingpool[kk,:] = Ft2

	return matingpool

def crossover_process(matingpool,pcross):
    npop,nvar = matingpool.shape
    population = np.zeros((npop,nvar))
    for jj in range(0,npop,2):
        idx1 = np.random.randint(npop)
        idx2 = np.random.randint(npop)
        if idx2 == idx1:
            while idx2 == idx1:
                idx2 = np.random.randint(npop)
        n = np.random.rand()
        p1 = matingpool[idx1-1,:]
        p2 = matingpool[idx2-1,:]
        if(n < pcross):
            child = SBX(p1,p2,nvar)
            population[jj,:] = child[0]
            population[jj+1,:] = child[1]
        else:
            population[jj,:] = p1
            population[jj+1,:] = p2

    return population

def polymut(var,nvar,pmut):
	nm = 20
	mutchrom = var
	const = 1/(1+nm)
	for i in range(nvar):
		u = np.random.rand()
		rn = np.random.rand()
		if(rn < pmut):
			delta = min(var[i],1-var[i])
			if(u <= 0.5):
				delta2 = (2*u+(1-2*u)*(1-delta)**(nm+1))**const - 1
			else:
				delta2 = 1 - (2*(1-u)+2*(u-0.5)*(1-delta)**(nm+1))**const
			mutchrom[i] = mutchrom[i] + delta2
		if (mutchrom[i] < 0):
			mutchrom[i] = 0
		if (mutchrom[i] > 1):
			mutchrom[i] = 1
	return mutchrom

def mutation_process(population,pmut):
	npop,nvar = population.shape
	population_mut = population.copy()
	for kk in range(npop):
		var = population[kk,:]
		population_mut[kk,:] = polymut(var,nvar,pmut)

	return population_mut

def calc_denorm(individual,lb,ub):
	indi_denorm = lb + (ub-lb)*individual
	return indi_denorm

def SBX(indi1,indi2,nvar):
    nc = 20
    c1 = np.zeros(nvar)
    c2 = np.zeros(nvar)
    count = 0
    for i in range(nvar):
        nr = np.random.rand()
        if nr < 0.5:
            u = np.random.rand()
            sn = 1e-16
            if abs(indi1[i]-indi2[i]) > sn:
                if indi1[i] < indi2[i]:
                    p1 = indi1[i]
                    p2 = indi2[i]
                else:
                    p1 = indi2[i]
                    p2 = indi1[i]

                const1 = p1
                const2 = 1-p2

                beta1 = 1 + 2*const1/(p2-p1)
                beta2 = 1 + 2*const2/(p2-p1)
                alpha1 = 2 - beta1**(-(nc+1))
                alpha2 = 2 - beta2**(-(nc+1))

                if u <= (1/alpha1):
                    beta_l1 = (alpha1*u)**(1/(nc+1))
                else:
                    beta_l1 = (1/(2-alpha1*u))**(1/(nc+1))

                c1[i] = 0.5*((p1+p2)-beta_l1*abs(p2-p1))

                if u <= (1/alpha2):
                    beta_l2 = (alpha2*u)**(1/(nc+1))
                else:
                    beta_l2 = (1/(2-alpha2*u))**(1/(nc+1))

                c2[i] = 0.5*((p1+p2)+beta_l2*abs(p2-p1))

            else:
                c1[i] = indi1[i]
                c2[i] = indi2[i]
        else:
            c1[i] = indi1[i]
            c2[i] = indi2[i]
        if(c1[i] < 0 or c1[i] > 1) or (c2[i] < 0 or c2[i] > 1):
            count = count+1
        if(count > 0):
            print(count)

    c = [c1,c2]
    return c