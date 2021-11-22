def EPF_quantile(data,var_name,probs = [i*0.01 for i in range(101)],J = [4,8,8],var_range = None):
  
    import numpy
    import pandas
    import math
    import time
    import scipy.optimize

    def data2var_scale(data_rdd,var_name):

        def f(x,var_name):
            res = []
            res.append(('-minimum',-x[var_name]))
            res.append(('maximum',x[var_name]))            
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,var_name)).reduceByKey(lambda x,y: max(x,y)).take(2)

        keys = list(map(lambda x:x[0],mr_out))
        values = list(map(lambda x:x[1],mr_out))

        for i in range(len(values)):
            if keys[i] == '-minimum':
                minimum = -values[i]
            else:
                maximum = values[i]

        M = 1/(maximum-minimum)
        C = -minimum/(maximum-minimum)

        var_scale = {'M' : M,'C' : C}        

        return(var_scale)

    def data2quantile_WEP_3_parameter(data_rdd,var_name,J_1,J_2,J_3,M,C):
        K = 2*J_1*J_2*J_3+1

        def get_quantile_summand_3_parameter(x,J_1,J_2,J_3):
            import math
            T = [(0,1)]
            c_1 = math.cos(x)
            c_2 = math.cos(2*J_1*x)
            c_3 = math.cos(2*J_1*J_2*x)
            m_1 = c_1*c_1
            c = c_1
            s_1 = math.sqrt(1-m_1)
            s = s_1
            m_2 = 1
            m_3 = 1
            I = 1

            for j_3 in range(J_3-1):
                for j_2 in range(J_2):
                    for j_1 in range(J_1):
                        T.append((I,c))
                        c *= m_1
                        I += 1
                        T.append((I,s))
                        s *= m_1
                        I += 1
                    m_2 *= c_2
                    c = c_1*m_2
                    s = s_1*m_2
                m_3 *= c_3
                m_2 = m_3
                c = c_1*m_2
                s = s_1*m_2

            for j_2 in range(J_2-1):
                for j_1 in range(J_1):
                    T.append((I,c))
                    c *= m_1
                    I += 1
                    T.append((I,s))
                    s *= m_1
                    I += 1 
                m_2 *= c_2
                c = c_1*m_2
                s = s_1*m_2

            for j_1 in range(J_1-1):
                T.append((I,c))
                c *= m_1
                I += 1
                T.append((I,s))
                s *= m_1
                I += 1

            T.append((I,c))
            I +=1
            T.append((I,s))  

            return(T)

        def f(x,J_1,J_2,J_3,var_name,M,C):
            x_scaled = x[var_name]*M+C
            res = get_quantile_summand_3_parameter(x_scaled,J_1,J_2,J_3)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,J_1,J_2,J_3,var_name,M,C)).reduceByKey(lambda x,y: x+y).take(K)

        Trig_3D_power_SEP = [None]*K
        Trig_3D_power_WEP = []

        for k in range(K):
            Trig_3D_power_SEP[mr_out[k][0]] = mr_out[k][1]
        for k in range(K):
            Trig_3D_power_WEP.append(Trig_3D_power_SEP[k]/Trig_3D_power_SEP[0])

        Trig_WEP = Trig_3D_power_WEP2Trig_WEP(Trig_3D_power_WEP,J_1,J_2,J_3)

        return(Trig_WEP)

    def data2quantile_WEP_4_parameter(data_rdd,var_name,J_1,J_2,J_3,J_4,M,C):
        K = 2*J_1*J_2*J_3*J_4+1

        def get_quantile_summand_4_parameter(x,J_1,J_2,J_3,J_4):
            T = [(0,1)]
            c_1 = math.cos(x)
            c_2 = math.cos(2*J_1*x)
            c_3 = math.cos(2*J_1*J_2*x)
            c_4 = math.cos(2*J_1*J_2*J_3*x)
            m_1 = c_1*c_1
            c = c_1
            s_1 = math.sqrt(1-m_1)
            s = s_1
            m_2 = 1
            m_3 = 1
            m_4 = 1
            I = 1

            for j_4 in range(J_4-1):
                for j_3 in range(J_3):
                    for j_2 in range(J_2):
                        for j_1 in range(J_1):
                            T.append((I,c))
                            c *= m_1
                            I += 1
                            T.append((I,s))
                            s *= m_1
                            I += 1
                        m_2 *= c_2
                        c = c_1*m_2
                        s = s_1*m_2
                    m_3 *= c_3
                    m_2 = m_3
                    c = c_1*m_2
                    s = s_1*m_2
                m_4 *= c_4
                m_3 = m_4
                m_2 = m_3
                c = c_1*m_2
                s = s_1*m_2

            for j_3 in range(J_3-1):
                for j_2 in range(J_2):
                    for j_1 in range(J_1):
                        T.append((I,c))
                        c *= m_1
                        I += 1
                        T.append((I,s))
                        s *= m_1
                        I += 1
                    m_2 *= c_2
                    c = c_1*m_2
                    s = s_1*m_2
                m_3 *= c_3
                m_2 = m_3
                c = c_1*m_2
                s = s_1*m_2

            for j_2 in range(J_2-1):
                for j_1 in range(J_1):
                    T.append((I,c))
                    c *= m_1
                    I += 1
                    T.append((I,s))
                    s *= m_1
                    I += 1 
                m_2 *= c_2
                c = c_1*m_2
                s = s_1*m_2

            for j_1 in range(J_1-1):
                T.append((I,c))
                c *= m_1
                I += 1
                T.append((I,s))
                s *= m_1
                I += 1

            T.append((I,c))
            I +=1
            T.append((I,s))  

            return(T)

        def f(x,J_1,J_2,J_3,J_4,var_name,M,C):
            x_scaled = x[var_name]*M+C
            res = get_quantile_summand_4_parameter(x_scaled,J_1,J_2,J_3,J_4)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,J_1,J_2,J_3,J_4,var_name,M,C)).reduceByKey(lambda x,y: x+y).take(K)

        Trig_4D_power_SEP = [None]*(K)
        Trig_4D_power_WEP = []

        for k in range(K):
            Trig_4D_power_SEP[mr_out[k][0]] = mr_out[k][1]
        for k in range(K):
            Trig_4D_power_WEP.append(Trig_4D_power_SEP[k]/Trig_4D_power_SEP[0])

        Trig_WEP = Trig_4D_power_WEP2Trig_WEP(Trig_4D_power_WEP,J_1,J_2,J_3,J_4)

        return(Trig_WEP)

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

    def Cosine_4D_product_WEP2Cosine_WEP(WEP):

        J = WEP.shape
        T = WEP.reshape((J[0],J[1]*J[2]*J[3]))
        T[0] = Cosine_3D_product_WEP2Cosine_WEP(WEP[0])
        add_T = T[0]

        for j in range(1,J[0]):
            add_T = 2*Cosine_3D_product_WEP2Cosine_WEP(WEP[j])-add_T[::-1]
            T[j] = add_T

        return(T.reshape((J[0]*J[1]*J[2]*J[3])))    

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

    def Sine_4D_product_WEP2Sine_WEP(WEP):
        J = WEP.shape
        T = WEP.reshape((J[0],J[1]*J[2]*J[3]))
        T[0] = Sine_3D_product_WEP2Sine_WEP(WEP[0])
        add_T = T[0]

        for j in range(1,J[0]):
            add_T = 2*Sine_3D_product_WEP2Sine_WEP(WEP[j])+add_T[::-1]
            T[j] = add_T

        return(T.reshape((J[0]*J[1]*J[2]*J[3])))    

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

    def Trig_4D_power_WEP2Trig_WEP(WEP,J_1,J_2,J_3,J_4):
        WEP_0 = WEP[0]
        WEP = numpy.delete(WEP,0)
        WEP = numpy.rollaxis(WEP.reshape((J_1*J_2*J_3*J_4,2)),1)

        Cosine_4D_power_WEP = WEP[0].reshape(J_4,J_3,J_2,J_1)
        Cosine_4D_temporary_WEP_1 = numpy.apply_along_axis(ChevTrA,3,Cosine_4D_power_WEP)
        Cosine_4D_temporary_WEP_2 = numpy.apply_along_axis(ChevTrB,2,Cosine_4D_temporary_WEP_1)
        Cosine_4D_temporary_WEP_3 = numpy.apply_along_axis(ChevTrB,1,Cosine_4D_temporary_WEP_2)
        Cosine_4D_product_WEP = numpy.apply_along_axis(ChevTrB,0,Cosine_4D_temporary_WEP_3)
        Cosine_WEP = Cosine_4D_product_WEP2Cosine_WEP(Cosine_4D_product_WEP)
        
        Sine_4D_power_WEP = WEP[1].reshape(J_4,J_3,J_2,J_1)
        Sine_4D_temporary_WEP_1 = numpy.apply_along_axis(SineTr,3,Sine_4D_power_WEP)
        Sine_4D_temporary_WEP_2 = numpy.apply_along_axis(ChevTrA,3,Sine_4D_temporary_WEP_1)
        Sine_4D_temporary_WEP_3 = numpy.apply_along_axis(OddNeg,3,Sine_4D_temporary_WEP_2)
        Sine_4D_temporary_WEP_4 = numpy.apply_along_axis(ChevTrB,2,Sine_4D_temporary_WEP_3)
        Sine_4D_temporary_WEP_5 = numpy.apply_along_axis(ChevTrB,1,Sine_4D_temporary_WEP_4)
        Sine_4D_product_WEP = numpy.apply_along_axis(ChevTrB,0,Sine_4D_temporary_WEP_5)
        Sine_WEP = Sine_4D_product_WEP2Sine_WEP(Sine_4D_product_WEP)

        Trig_WEP = numpy.append(WEP_0,numpy.rollaxis(numpy.vstack([Cosine_WEP,Sine_WEP]),1).flatten())

        return(Trig_WEP)

    def get_quantile_multiplier(m,J):
        F = numpy.zeros(2*J+1)
        m_inc = 2*m
        F[0] = .5

        for j in range(1,2*J,2):
            F[j] = (2/(math.pi*j))*math.sin(m)
            F[j+1] = (-2/(math.pi*j))*math.cos(m)
            m += m_inc

        return(F)

    def WEP2quantile(WEP,probs,M,C):
        J = int((len(WEP)-1)/2)
        Q = [None]*len(probs)

        for i in range(len(probs)):
            if probs[i] == 0:
                Q[i] = -C/M
            elif probs[i] == 1:
                Q[i] = (1-C)/M
            else:
                Q[i] = (scipy.optimize.brentq(lambda m: WEP.dot(get_quantile_multiplier(m,J))-probs[i],0,1)-C)/M

        return(Q)
    
    start_time = time.time()
        
    if (not isinstance(data,pyspark.sql.dataframe.DataFrame)):
        raise ValueError('"data" should be a pyspark data-frame object!')
    
    if (not var_name in data.columns) or (data.dtypes[data.columns.index(var_name)][1] != 'double'):
        raise ValueError('"var_name" should be a character representing a numeric-variable present in the data!')
        
    if (not isinstance(probs,list)) or (max(probs) > 1) or (min(probs) < 0):
        raise ValueError('"probs" should be a list of f-values!')
        
    if (var_range != None):
        if (isinstance(var_range,list) == False) or (len(var_range) != 2) or (var_range[0] > var_range[1]):
            raise ValueError('"var_range" should be a list containing the minimum and the maximum of the numeric-variable!')
        else:
            M = 1/(var_range[1]-var_range[0])
            C = -var_range[0]/(var_range[1]-var_range[0])
    else:
        scales = data2var_scale(data.rdd,var_name)
        M = scales['M']
        C = scales['C']
        
    if (not isinstance(J,list)) or (not all(isinstance(j, int) for j in J)) or (not len(J) in [3,4]):
        raise ValueError('"J" should be a list containing 3 or 4 integer parameters!')
    
    if (len(J) == 3):
        WEP = data2quantile_WEP_3_parameter(data.rdd,var_name,J[0],J[1],J[2],M,C)
    
    if (len(J) == 4):
        WEP = data2quantile_WEP_4_parameter(data.rdd,var_name,J[0],J[1],J[2],J[3],M,C)
                                                                                                 
    Q = WEP2quantile(WEP,probs,M,C)
    
    stop_time = time.time()
        
    EPF_quantile_output = {'Quantiles' : pandas.DataFrame({'f_values' : probs,'Quantile_values' : Q}),'Computation_time( in seconds)' : stop_time-start_time}
    
    return(EPF_quantile_output)
