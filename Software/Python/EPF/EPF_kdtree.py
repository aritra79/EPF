def EPF_kdtree(data, var_names, D = 3, J = [3,6]):
    
    import numpy
    import pandas
    import math
    import time
    import scipy.optimize

    def data2var_scales(data_rdd, var_names):
        P = len(var_names)

        def f(x,var_names):
            res = []
            for i in range(P):
                res.append(((i,'-minimum'),-x[var_names[i]]))
                res.append(((i,'maximum'),x[var_names[i]]))            
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,var_names)).reduceByKey(lambda x,y: max(x,y)).take(2*len(var_names))

        keys_var_id = list(map(lambda x:x[0][0],mr_out))
        keys_type = list(map(lambda x:x[0][1],mr_out))
        values = list(map(lambda x:x[1],mr_out))
        minimum = [None]*P
        maximum = [None]*P

        for i in range(len(values)):
            if keys_type[i] == '-minimum':
                minimum[keys_var_id[i]] = -values[i]
            else:
                maximum[keys_var_id[i]] = values[i]

        M = []
        C = []

        for i in range(P):
            M.append(1/(maximum[i]-minimum[i]))
            C.append(-minimum[i]/(maximum[i]-minimum[i]))

        var_scales = {'M' : M, 'C' : C}        

        return(var_scales)

    def data2kdtree_WEP_2_parameter(data_rdd, var_names, J_1, J_2, M, C):
        import numpy
        P = len(var_names)
        J = J_1*J_2
        K = (2*J+1)**P

        def get_kdtree_summand_2_parameter(x, J_1, J_2, J, P):
            import math
            T = []
            M = []
            C = []

            for p in range(P-1):
                C.append(1)
                c_1 = math.cos(x[p])
                c_2 = math.cos(2*J_1*x[p])
                m_1 = c_1*c_1
                c = c_1
                s_1 = math.sqrt(1-m_1)
                s = s_1
                m_2 = 1

                for j_2 in range(J_2-1):
                    for j_1 in range(J_1):
                        C.append(c)
                        c *= m_1
                        C.append(s)
                        s *= m_1 
                    m_2 *= c_2
                    c = c_1*m_2
                    s = s_1*m_2

                for j_1 in range(J_1-1):
                    C.append(c)
                    c *= m_1
                    C.append(s)
                    s *= m_1

                C.append(c)
                C.append(s)  


            C.append(1)
            c_1 = math.cos(x[P-1])
            c_2 = math.cos(2*J_1*x[P-1])
            m_1 = c_1*c_1
            c = c_1
            s_1 = math.sqrt(1-m_1)
            s = s_1
            m_2 = 1

            for j_2 in range(J_2-1):
                for j_1 in range(J_1):
                    C.append(c)
                    c *= m_1
                    C.append(s)
                    s *= m_1 
                m_2 *= c_2
                c = c_1*m_2
                s = s_1*m_2

            for j_1 in range(J_1-1):
                C.append(c)
                c *= m_1
                C.append(s)
                s *= m_1

            C.append(c)
            C.append(s)

            for c in range(2*J+1):
                M.append(C[c])
            K = c+1

            if P > 2:
                for p in range(1,P-1):
                    T = []
                    for j in range(2*J+1):
                        c += 1 
                        for k in range(K):
                            T.append(M[k]*C[c])       
                    K *= 2*J+1
                    M = []
                    for k in range(K):
                        M.append(T[k])

            T = []
            t = 0
            for j in range(2*J+1):
                c += 1
                for k in range(K):
                    T.append((t,M[k]*C[c]))
                    t += 1

            return(T)    

        def f(x, var_names, M, C, J_1, J_2, J, P):
            P = len(var_names)
            x_scaled = []
            for i in range(P):
                x_scaled.append(x[var_names[i]]*M[i]+C[i])
            res = get_kdtree_summand_2_parameter(x_scaled, J_1, J_2, J, P)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x, var_names, M, C, J_1, J_2, J, P)).reduceByKey(lambda x,y: x+y).take(K)

        kdtree_SEP = numpy.zeros(K)
        kdtree_WEP = numpy.zeros(K)

        for k in range(K):
            kdtree_SEP[mr_out[k][0]] = mr_out[k][1]
        for k in range(K):
            kdtree_WEP[k] = kdtree_SEP[k]/kdtree_SEP[0]

        kdtree_WEP = kdtree_WEP.reshape(numpy.repeat(2*J+1,P))

        def f_tr(WEP):
            res = Trig_2D_power_WEP2Trig_WEP(WEP, J_1, J_2)
            return(res)
        
        for p in range(P):
            kdtree_WEP = numpy.apply_along_axis(f_tr, p, kdtree_WEP)

        return(kdtree_WEP)

    def data2kdtree_WEP_3_parameter(data_rdd, var_names, J_1, J_2, J_3, M, C):
        import numpy
        P = len(var_names)
        J = J_1*J_2*J_3
        K = (2*J+1)**P

        def get_kdtree_summand_3_parameter(x, J_1, J_2, J_3, J, P):
            import math
            T = []
            M = []
            C = []

            for p in range(P-1):
                C.append(1)
                c_1 = math.cos(x[p])
                c_2 = math.cos(2*J_1*x[p])
                c_3 = math.cos(2*J_1*J_2*x[p])
                m_1 = c_1*c_1
                c = c_1
                s_1 = math.sqrt(1-m_1)
                s = s_1
                m_2 = 1
                m_3 = 1

                for j_3 in range(J_3-1):
                    for j_2 in range(J_2):
                        for j_1 in range(J_1):
                            C.append(c)
                            c *= m_1
                            C.append(s)
                            s *= m_1
                        m_2 *= c_2
                        c = c_1*m_2
                        s = s_1*m_2
                    m_3 *= c_3
                    m_2 = m_3
                    c = c_1*m_2
                    s = s_1*m_2

                for j_2 in range(J_2-1):
                    for j_1 in range(J_1):
                        C.append(c)
                        c *= m_1
                        C.append(s)
                        s *= m_1 
                    m_2 *= c_2
                    c = c_1*m_2
                    s = s_1*m_2

                for j_1 in range(J_1-1):
                    C.append(c)
                    c *= m_1
                    C.append(s)
                    s *= m_1

                C.append(c)
                C.append(s)  


            C.append(1)
            c_1 = math.cos(x[P-1])
            c_2 = math.cos(2*J_1*x[P-1])
            c_3 = math.cos(2*J_1*J_2*x[P-1])
            m_1 = c_1*c_1
            c = c_1
            s_1 = math.sqrt(1-m_1)
            s = s_1
            m_2 = 1
            m_3 = 1

            for j_3 in range(J_3-1):
                for j_2 in range(J_2):
                    for j_1 in range(J_1):
                        C.append(c)
                        c *= m_1
                        C.append(s)
                        s *= m_1
                    m_2 *= c_2
                    c = c_1*m_2
                    s = s_1*m_2
                m_3 *= c_3
                m_2 = m_3
                c = c_1*m_2
                s = s_1*m_2

            for j_2 in range(J_2-1):
                for j_1 in range(J_1):
                    C.append(c)
                    c *= m_1
                    C.append(s)
                    s *= m_1 
                m_2 *= c_2
                c = c_1*m_2
                s = s_1*m_2

            for j_1 in range(J_1-1):
                C.append(c)
                c *= m_1
                C.append(s)
                s *= m_1

            C.append(c)
            C.append(s)

            for c in range(2*J+1):
                M.append(C[c])
            K = c+1

            if P > 2:
                for p in range(1,P-1):
                    T = []
                    for j in range(2*J+1):
                        c += 1 
                        for k in range(K):
                            T.append(M[k]*C[c])       
                    K *= 2*J+1
                    M = []
                    for k in range(K):
                        M.append(T[k])

            T = []
            t = 0
            for j in range(2*J+1):
                c += 1
                for k in range(K):
                    T.append((t,M[k]*C[c]))
                    t += 1

            return(T)    

        def f(x, var_names, M, C, J_1, J_2, J_3, J, P):
            P = len(var_names)
            x_scaled = []
            for i in range(P):
                x_scaled.append(x[var_names[i]]*M[i]+C[i])
            res = get_kdtree_summand_3_parameter(x_scaled, J_1, J_2, J_3, J, P)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x, var_names, M, C, J_1, J_2, J_3, J, P)).reduceByKey(lambda x,y: x+y).take(K)

        kdtree_SEP = numpy.zeros(K)
        kdtree_WEP = numpy.zeros(K)

        for k in range(K):
            kdtree_SEP[mr_out[k][0]] = mr_out[k][1]
        for k in range(K):
            kdtree_WEP[k] = kdtree_SEP[k]/kdtree_SEP[0]

        kdtree_WEP = kdtree_WEP.reshape(numpy.repeat(2*J+1,P))

        def f_tr(WEP):
            res = Trig_3D_power_WEP2Trig_WEP(WEP, J_1, J_2, J_3)
            return(res)

        for p in range(P):
            kdtree_WEP = numpy.apply_along_axis(f_tr, p, kdtree_WEP)

        return(kdtree_WEP)
    
    def ChevTrA(T):

        J = len(T)
        E = numpy.zeros(J)
        O = numpy.zeros(J)

        for j in range(J):
            E[j] = T[j]
            O[j] = T[j]

        for j in range(J):
            T[j] = O[0]

            for k in range(J-j-1):
                E[k] = 2*O[k+1]-E[k]
                O[k] = 2*E[k]-O[k]

        return(T)

    def ChevTrB(T):

        J = len(T)
        E = numpy.zeros(J)
        O = numpy.zeros(J-1)
        E[0] = T[0]

        for j in range(1,J):
            E[j] = T[j]
            O[j-1] = T[j] 

        for j in range(math.floor(J/2)):
            T[2*j] = E[0]
            T[2*j+1] = O[0]

            for k in range(J-2*j-2):
                E[k] = 2*O[k+1]-E[k]

            for k in range(J-2*j-3):
                O[k] = 2*E[k+1]-O[k]

        if J % 2 == 1:
            T[J-1] = E[0]   

        return(T)

    def SineTr(T):

        J = len(T)
        S = numpy.zeros(J)

        for k in range(1,J):
            for j in range(k,J):
                S[j] = T[j-1]-T[j]

            for j in range(k,J):
                T[j] = S[j]

        return(T)

    def OddNeg(T):

        for j in range(len(T)):
            T[j] = (-1)**j*T[j]

        return(T)

    def Cosine_2D_product_WEP2Cosine_WEP(WEP):

        J = WEP.shape
        T = WEP
        add_T = T[0]

        for j in range(1,J[0]):
            add_T = 2*WEP[j]-add_T[::-1]
            T[j] = add_T

        return(T.reshape((J[0]*J[1])))

    def Cosine_3D_product_WEP2Cosine_WEP(WEP):

        J = WEP.shape
        T = WEP.reshape((J[0],J[1]*J[2]))
        T[0] = Cosine_2D_product_WEP2Cosine_WEP(WEP[0])
        add_T = T[0]

        for j in range(1,J[0]):
            add_T = 2*Cosine_2D_product_WEP2Cosine_WEP(WEP[j])-add_T[::-1]
            T[j] = add_T

        return(T.reshape((J[0]*J[1]*J[2])))

    def Sine_2D_product_WEP2Sine_WEP(WEP):
        J = WEP.shape
        T = WEP
        add_T = T[0]

        for j in range(1,J[0]):
            add_T = 2*WEP[j]+add_T[::-1]
            T[j] = add_T

        return(T.reshape((J[0]*J[1])))

    def Sine_3D_product_WEP2Sine_WEP(WEP):
        J = WEP.shape
        T = WEP.reshape((J[0],J[1]*J[2]))
        T[0] = Sine_2D_product_WEP2Sine_WEP(WEP[0])
        add_T = T[0]

        for j in range(1,J[0]):
            add_T = 2*Sine_2D_product_WEP2Sine_WEP(WEP[j])+add_T[::-1]
            T[j] = add_T

        return(T.reshape((J[0]*J[1]*J[2])))    

    def Trig_2D_power_WEP2Trig_WEP(WEP, J_1, J_2):
        WEP_0 = WEP[0]
        WEP = numpy.delete(WEP,0)
        WEP = numpy.rollaxis(WEP.reshape((J_1*J_2,2)),1)

        Cosine_2D_power_WEP = WEP[0].reshape(J_2,J_1)
        Cosine_2D_temporary_WEP = numpy.apply_along_axis(ChevTrA,1,Cosine_2D_power_WEP)
        Cosine_2D_product_WEP = numpy.apply_along_axis(ChevTrB,0,Cosine_2D_temporary_WEP)
        Cosine_WEP = Cosine_2D_product_WEP2Cosine_WEP(Cosine_2D_product_WEP)
        
        Sine_2D_power_WEP = WEP[1].reshape(J_2,J_1)
        Sine_2D_temporary_WEP_1 = numpy.apply_along_axis(SineTr,1,Sine_2D_power_WEP)
        Sine_2D_temporary_WEP_2 = numpy.apply_along_axis(ChevTrA,1,Sine_2D_temporary_WEP_1)
        Sine_2D_temporary_WEP_3 = numpy.apply_along_axis(OddNeg,1,Sine_2D_temporary_WEP_2)
        Sine_2D_product_WEP = numpy.apply_along_axis(ChevTrB,0,Sine_2D_temporary_WEP_3)
        Sine_WEP = Sine_2D_product_WEP2Sine_WEP(Sine_2D_product_WEP)
        
        Trig_WEP = numpy.append(WEP_0,numpy.rollaxis(numpy.vstack([Cosine_WEP,Sine_WEP]),1).flatten())

        return(Trig_WEP)    
    
    def Trig_3D_power_WEP2Trig_WEP(WEP,J_1,J_2,J_3):
        WEP_0 = WEP[0]
        WEP = numpy.delete(WEP,0)
        WEP = numpy.rollaxis(WEP.reshape((J_1*J_2*J_3,2)),1)

        Cosine_3D_power_WEP = WEP[0].reshape(J_3,J_2,J_1)
        Cosine_3D_temporary_WEP_1 = numpy.apply_along_axis(ChevTrA,2,Cosine_3D_power_WEP)
        Cosine_3D_temporary_WEP_2 = numpy.apply_along_axis(ChevTrB,1,Cosine_3D_temporary_WEP_1)
        Cosine_3D_product_WEP = numpy.apply_along_axis(ChevTrB,0,Cosine_3D_temporary_WEP_2)
        Cosine_WEP = Cosine_3D_product_WEP2Cosine_WEP(Cosine_3D_product_WEP)

        Sine_3D_power_WEP = WEP[1].reshape(J_3,J_2,J_1)
        Sine_3D_temporary_WEP_1 = numpy.apply_along_axis(SineTr,2,Sine_3D_power_WEP)
        Sine_3D_temporary_WEP_2 = numpy.apply_along_axis(ChevTrA,2,Sine_3D_temporary_WEP_1)
        Sine_3D_temporary_WEP_3 = numpy.apply_along_axis(OddNeg,2,Sine_3D_temporary_WEP_2)
        Sine_3D_temporary_WEP_4 = numpy.apply_along_axis(ChevTrB,1,Sine_3D_temporary_WEP_3)
        Sine_3D_product_WEP = numpy.apply_along_axis(ChevTrB,0,Sine_3D_temporary_WEP_4)
        Sine_WEP = Sine_3D_product_WEP2Sine_WEP(Sine_3D_product_WEP)

        Trig_WEP = numpy.append(WEP_0,numpy.rollaxis(numpy.vstack([Cosine_WEP,Sine_WEP]),1).flatten())

        return(Trig_WEP)    

    def get_kdtree_multiplier(a, b, J):
        F = numpy.zeros(2*J+1)
        a_inc = 2*a
        b_inc = 2*b
        F[0] = 1

        if a > 0:
            F[0] -= .5
            a_inc = 2*a

            for j in range(1,2*J,2):
                F[j] -= 2/(math.pi*j)*math.sin(a)
                F[j+1] += 2/(math.pi*j)*math.cos(a)
                a += a_inc

        if b < 1:
            F[0] -= .5
            a_inc = 2*a

            for j in range(1,2*J,2):
                F[j] += 2/(math.pi*j)*math.sin(b)
                F[j+1] -= 2/(math.pi*j)*math.cos(b)
                b += b_inc

        return(F)

    def WEP2multiplier_info(WEP, lower, upper, d):
        P = len(WEP.shape)
        d = d%P
        J = int((WEP.shape[1]-1)/2)
        multiplier_array = 1

        for p in range(P):
            if p != d:
                multiplier_array = numpy.outer(multiplier_array, get_kdtree_multiplier(lower[p], upper[p], J))

        multiplier_array = multiplier_array.flatten()
        multiplier = numpy.rollaxis(WEP,P-d-1).reshape(2*J+1,(2*J+1)**(P-1)).dot(multiplier_array)
        multiplier_info = {'d' : d, 'multiplier' : multiplier}

        return(multiplier_info)

    def WEP2median_info(WEP, lower, upper, d):
        multiplier_info = WEP2multiplier_info(WEP, lower, upper, d)
        d = multiplier_info['d']
        multiplier = multiplier_info['multiplier']
        J = int((WEP.shape[1]-1)/2)

        def fn(m):
            return(multiplier.dot(get_kdtree_multiplier(lower[d],m,J)-get_kdtree_multiplier(m,upper[d],J)))

        median = scipy.optimize.brentq(lambda m: fn(m),lower[d],upper[d])
        median_info = {'d' : d, 'median' : median}

        return(median_info)

    def WEP2sub_neighborhoods_unscaled_info(WEP, neighborhoods_unscaled, d, M, C):
        P = len(WEP.shape)
        J = int((WEP.shape[1]-1)/2)
        I = neighborhoods_unscaled.shape[0]
        new_neighborhoods_unscaled = numpy.zeros((2*I,2,P))
        new_medians = numpy.zeros(I)

        for i in range(I):
            median_info = WEP2median_info(WEP = WEP, lower = neighborhoods_unscaled[i][0], upper = neighborhoods_unscaled[i][1], d = d)
            d = median_info['d']
            new_medians[i] = (median_info['median']-C[d])/M[d]
            new_neighborhoods_unscaled[2*i][0] = neighborhoods_unscaled[i][0]
            new_neighborhoods_unscaled[2*i][1] = neighborhoods_unscaled[i][1]
            new_neighborhoods_unscaled[2*i+1][0] = neighborhoods_unscaled[i][0]
            new_neighborhoods_unscaled[2*i+1][1] = neighborhoods_unscaled[i][1]
            new_neighborhoods_unscaled[2*i][1][d] = median_info['median']
            new_neighborhoods_unscaled[2*i+1][0][d] = median_info['median']

        neighborhoods_unscaled_info = {'neighborhoods_unscaled' : new_neighborhoods_unscaled, 'medians' : new_medians}

        return(neighborhoods_unscaled_info)

    def scale_neighborhoods(neighborhoods_unscaled, M, C):
        [I,J,P] = neighborhoods_unscaled.shape
        neighborhoods_scaled = numpy.zeros((I,J,P))

        for i in range(I):
            for j in range(J):
                for p in range(P):
                    neighborhoods_scaled[i][j][p] = (neighborhoods_unscaled[i][j][p]-C[p])/M[p]            

        return(neighborhoods_scaled)        

    def WEP2kdtree_info(WEP, D, M, C):
        P = len(WEP.shape)
        neighborhoods_unscaled = numpy.zeros((1,2,P))
        neighborhoods_unscaled[0][1] = numpy.ones(P)
        medians = [None]*D

        for d in range(D):
            neighborhoods_unscaled_info = WEP2sub_neighborhoods_unscaled_info(WEP = WEP, neighborhoods_unscaled = neighborhoods_unscaled, d = d, M = M, C = C)
            neighborhoods_unscaled = neighborhoods_unscaled_info['neighborhoods_unscaled']
            medians[d] = neighborhoods_unscaled_info['medians']

        neighborhoods = scale_neighborhoods(neighborhoods_unscaled = neighborhoods_unscaled, M = M, C = C)
        kdtree_info = {'neighborhoods' : neighborhoods, 'medians' : medians}

        return(kdtree_info)

    def kdtree_medians2cellcounts(data_rdd, var_names, medians):
        K = 2**len(medians)

        def f(x, var_names, medians):
            P = len(var_names)
            D = len(medians)
            i = 0
            for d in range(D):
                dd = d%P
                if x[var_names[dd]] == medians[d][i]:
                    i = -1
                    break
                elif x[var_names[dd]] < medians[d][i]:
                    i = 2*i
                else:
                    i = 2*i+1
            res = (i,1)        
            return(res)    

        mr_out = data_rdd.map(lambda x: f(x,var_names, medians)).reduceByKey(lambda x,y: x+y).take(K+1)

        cellcounts = numpy.zeros(K)

        for k in range(K):
            cellcounts[mr_out[k][0]] = mr_out[k][1]

        return(cellcounts)

    
    start_time = time.time()
        
    if (not isinstance(data, pyspark.sql.dataframe.DataFrame)):
        raise ValueError('"data" should be a pyspark data-frame object!')
    
    for var in var_names:
        if (not var in data.columns) or (data.dtypes[data.columns.index(var)][1] != 'double'):
            raise ValueError('Each element in "var_names" should be a character representing the explanatory-variable present in the data!')
    
    var_scales = data2var_scales(data.rdd, var_names)
        
    if ((type(D) != int) or (D < 2)):
        raise ValueError('"D" should be a positive integer representing the depth of the kdtree!')

    if (not isinstance(J,list)) or (not all(isinstance(j, int) for j in J)) or (not len(J) in [2,3]):
        raise ValueError('"J" should be a list containing 2 or 3 integer parameters!')
    
    if (len(J) == 2):
        WEP = data2kdtree_WEP_2_parameter(data.rdd, var_name, J[0], J[1], var_scales['M'], var_scales['C'])
    
    if (len(J) == 3):
        WEP = data2kdtree_WEP_3_parameter(data.rdd, var_name, J[0], J[1], J[2], var_scales['M'], var_scales['C'])
                                                                                                 
    kdtree_info = WEP2kdtree_info(WEP, D, var_scales['M'], var_scales['C'])
    tree = []
    
    for d in range(D):
        tree.append({'Level' : d, 'Splitting axis' : var_names[d%len(var_names)], 'Medians' : kdtree_info['medians'][d].tolist()})
        
    neighborhoods = []
    
    for i in range(kdtree_info['neighborhoods'].shape[0]):
        neighborhoods.append({'lower' : kdtree_info['neighborhoods'][i][0].tolist(), 'upper' : kdtree_info['neighborhoods'][i][1].tolist()})
    
    counts = kdtree_medians2cellcounts(data.rdd, var_names, kdtree_info['medians'])
    
    stop_time = time.time()
    
    kdtree = {'Tree' : tree, 'Neighborhoods' : neighborhoods, 'Counts' : counts, 'Computation_time_in_seconds' : stop_time-start_time}
    
    return(kdtree)

