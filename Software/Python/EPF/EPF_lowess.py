def EPF_lowess(data,var_x,var_y,x_0,J = [4,8,8],var_x_range = None,alpha = 0.01,K = 2):
    
    import numpy
    import pandas
    import math
    import time
    import scipy.optimize

    def data2var_scale(data_rdd,var_x):

        def f(x,var_x):
            res = []
            res.append(('-minimum',-x[var_x]))
            res.append(('maximum',x[var_x]))            
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,var_x)).reduceByKey(lambda x,y: max(x,y)).take(2)

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

    def data2lowess_WEP_3_parameter(data_rdd,var_x,J_1,J_2,J_3,M,C):
        K = 2*J_1*J_2*J_3+1

        def get_lowess_summand_3_parameter(x,J_1,J_2,J_3):
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

        def f(x,J_1,J_2,J_3,var_x,M,C):
            x_scaled = x[var_x]*M+C
            res = get_lowess_summand_3_parameter(x_scaled,J_1,J_2,J_3)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,J_1,J_2,J_3,var_x,M,C)).reduceByKey(lambda x,y: x+y).take(K)

        Trig_3D_power_SEP = [None]*K
        Trig_3D_power_WEP = []

        for k in range(K):
            Trig_3D_power_SEP[mr_out[k][0]] = mr_out[k][1]
        for k in range(K):
            Trig_3D_power_WEP.append(Trig_3D_power_SEP[k]/Trig_3D_power_SEP[0])

        Trig_WEP = Trig_3D_power_WEP2Trig_WEP(Trig_3D_power_WEP,J_1,J_2,J_3)

        return(Trig_WEP)

    def data2lowess_WEP_4_parameter(data_rdd,var_x,J_1,J_2,J_3,J_4,M,C):
        K = 2*J_1*J_2*J_3*J_4+1

        def get_lowess_summand_4_parameter(x,J_1,J_2,J_3,J_4):
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

        def f(x,J_1,J_2,J_3,J_4,var_x,M,C):
            x_scaled = x[var_x]*M+C
            res = get_lowess_summand_4_parameter(x_scaled,J_1,J_2,J_3,J_4)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,J_1,J_2,J_3,J_4,var_x,M,C)).reduceByKey(lambda x,y: x+y).take(K)

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

    def get_lowess_multiplier(x_0,h,J):
        F = numpy.zeros(2*J+1)
        x_0_inc = 2*x_0
        h_inc = 2*h       
        F[0] = 0

        for j in range(1,2*J,2):
            F[j] = (4/(math.pi*j))*math.cos(x_0)*math.sin(h)
            F[j+1] = (4/(math.pi*j))*math.sin(x_0)*math.sin(h)
            x_0 += x_0_inc
            h += h_inc

        return(F)

    def Trig_WEP2lowess_widths(WEP,x_0,alpha,M,C):
        J = int((len(WEP)-1)/2)
        h = [None]*len(x_0)

        for i in range(len(x_0)):
            xx = M*x_0[i]+C
            h[i] = scipy.optimize.brentq(lambda h: WEP.dot(get_lowess_multiplier(xx,h,J))-alpha,0,1)/M

        return(h)

    def get_lowess_powers(x,y,x_0,h,K):
        P = []
        I = 0
        for j in range(len(x_0)):
            if x_0[j] - h[j] < x and x_0[j] + h[j] > x:
                p = (1-abs((x-x_0[j])/h[j])**3)**3
                P.append((I,p))
                I += 1
                for k in range(1,2*K+1):
                    p *= x-x_0[j]
                    P.append((I,p))
                    I += 1
                p = y*(1-abs((x-x_0[j])/h[j])**3)**3
                P.append((I,p))
                I += 1
                for k in range(2*K+2,3*K+2):
                    p *= x-x_0[j]
                    P.append((I,p))
                    I += 1
            else:
                for k in range(3*K+2):
                    P.append((I,0))
                    I += 1
        return(P)  

    def data2lowess_predicted_values(data_rdd,var_x,var_y,x_0,h,K):

        def f(x,var_x,var_y,x_0,h,K):
            res = get_lowess_powers(x[var_x],x[var_y],x_0,h,K)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,var_x,var_y,x_0,h,K)).reduceByKey(lambda x,y: x+y).take(len(x_0)*(3*K+2))

        lowess_powers = numpy.zeros(len(x_0)*(3*K+2))

        for k in range(len(x_0)*(3*K+2)):
            lowess_powers[mr_out[k][0]] = mr_out[k][1]

        lowess_power_matrix = lowess_powers.reshape((len(x_0),3*K+2))
        y_0_hat = []

        for i in range(len(x_0)):
            y_0_hat.append(lowess_powers2predicted_value(lowess_power_matrix[i],K))

        return(y_0_hat)

    def lowess_powers2predicted_value(p):
        K = int((len(p)-2)/3)
        xx_power = numpy.zeros((K+1,K+1))

        for k in range(K+1):
            for l in range(K+1):
                xx_power[k,l] = p[k+l]

        xy_power = p[(2*K+1):(3*K+2)]
        y_0_hat = numpy.linalg.inv(xx_power).dot(xy_power)[0]

        return(y_0_hat)

    def get_lowess_match(x,y,x_0):
        T = []
        for i in range(len(x_0)):
            t = (i,[])
            if x_0[i] == x:
                t[1].append(y)
            T.append(t)
        return(T)             

    def data2residuals(data_rdd,var_x,var_y,x_0,y_0_hat):

        def f(x,var_x,var_y,x_0):
            res = get_lowess_match(x[var_x],x[var_y],x_0)
            return(res)

        mr_out = data_rdd.flatMap(lambda x: f(x,var_x,var_y,x_0)).reduceByKey(lambda x,y: x+y).take(len(x_0))

        y_actual = [None]*len(x_0)

        for i in range(len(x_0)):
            y_actual[mr_out[i][0]] = numpy.array(mr_out[i][1])

        lowess_stats = [None]*len(x_0)

        for i in range(len(x_0)):
            if len(y_actual[i]) > 0:
                lowess_stats[i] = {'x' : x_0[i],'y_actual_value' : y_actual[i].tolist(),'residual' : (y_actual[i]-y_0_hat[i]).tolist()}
            else:
                lowess_stats[i] = {'x' : x_0[i],'y_actual_value' : None,'residual' : None}

        return(lowess_stats)
        
    start_time = time.time()
        
    if (not isinstance(data,pyspark.sql.dataframe.DataFrame)):
        raise ValueError('"data" should be a pyspark data-frame object!')
    
    if (not var_x in data.columns) or (data.dtypes[data.columns.index(var_x)][1] != 'double'):
        raise ValueError('"var_x" should be a character representing the explanatory-variable present in the data!')
        
    if (not var_y in data.columns) or (data.dtypes[data.columns.index(var_y)][1] != 'double'):
        raise ValueError('"var_y" should be a character representing the response-variable present in the data!')
    
    if (var_x_range != None):
        if (isinstance(var_x_range,list) == False) or (len(var_x_range) != 2) or (var_x_range[0] > var_x_range[1]):
            raise ValueError('"var_x_range" should be a list containing the minimum and the maximum of the explanatory-variable!')
        else:
            M = 1/(var_x_range[1]-var_x_range[0])
            C = -var_x_range[0]/(var_x_range[1]-var_x_range[0])
    else:
        scales = data2var_scale(data.rdd,var_x)
        M = scales['M']
        C = scales['C']
        var_x_range = [-C/M,(1-C)/M]
        
    if (not isinstance(J,list)) or (not all(isinstance(j,int) for j in J)) or (not len(J) in [3,4]):
        raise ValueError('"J" should be a list containing 3 or 4 integer parameters!')
    
    if (len(J) == 3):
        WEP = data2lowess_WEP_3_parameter(data.rdd,var_x,J[0],J[1],J[2],M,C)
    
    if (len(J) == 4):
        WEP = data2lowess_WEP_4_parameter(data.rdd,var_x,J[0],J[1],J[2],J[3],M,C)
    
    if ((isinstance(x_0,list) == False) or (not all(isinstance(x,float) for x in x_0)) or (min(x_0) < var_x_range[0]) or (max(x_0) > var_x_range[1])):
        raise ValueError('"x_0" should be a list of test-data and it should have the same range as var_x!')
        
    if ((type(alpha) != float) or (alpha < 0) or (alpha > 1)):
        raise ValueError('"alpha" should be a numaber in (0,1)!')
        
    if ((type(K) != int) or (K < 1)):
        raise ValueError('"K" should be a positive integer!')
    
    h = Trig_WEP2lowess_width(WEP,scales['M'],scales['C'],x_0,alpha)
    y_hat = data2lowess_predicted_values(data.rdd,var_x,var_y,x_0,h,K)
    lowess_statistics = data2residuals(data.rdd,var_x,var_y,x_0,y_hat)
    error = numpy.array(sum(list(filter(lambda x : x != None,list(map(lambda x : x['residual'],lowess_statistics)))),[]))
    RSS = error.dot(error)
    
    stop_time = time.time()
    
    EPF_lowess_output = {'Lowess_fit' : pandas.DataFrame({'x' : x_0,'y_hat' : y_hat}),'Lowess_statistics' : lowess_statistics,'Residual_sum_of_squares' : RSS,'Computation_time_in_seconds' : stop_time-start_time}
    
    return(EPF_lowess_output)    
