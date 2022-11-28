import pyromat as pm
import numpy as np
import pylab as plt
import pandas as pd
import seaborn as sns


H2O = pm.get('mp.H2O')

plt.style.use('seaborn')

p1=8
p2 = 0.008
x1=1
nu_t_real=0.85
nu_p_real=0.85
Wc = 100e3


pm.config['unit_pressure'] = 'MPa'
pm.config['unit_temperature'] = 'C'


#CICLO RANKINE IDEAL

def Rankine(x=pm.get('mp.H2O'),p1=p1,p2=p2,x1=x1,x3=0,nu_p_real=nu_p_real,nu_t_real=nu_t_real,Wc = Wc): #Dejo por efecto los valores dados en el enunciado

        #Estado 1
        h1 = x.h(p=p1,x=x1)
        s1 = x.s(p=p1,x=x1)
        t1 = x.T(p=p1,x=x1)

        #Estado 2
        x2 = x.T_s(s=s1,p=p2, quality=True)[1]
        h2 = x.h(p=p2,x=x2)
        s2 = x.s(p=p2,x=x2)
        t2 = x.T(p=p2,x=x2)

        #Estado 2 real
        h2real = h1 - nu_t_real*(h1-h2)

        T2_real,x2_real  = H2O.T_h(h=h2real,p=p2, quality=True)
        s2_real  = H2O.s(T=T2_real,p=p2,x=x2_real)

        #Estado 3
        t3=t2
        x3=x3

        p3 = x.p(T=t3, x=x3)
        s3 = x.s(p=p3, T=t3, x=x3)
        h3 = x.h(p=p3, T=t3, x=x3)

        #Estado 4
        s4=s3
        p4=p1

        t4,x4 = x.T_s(s=s4,p=p4, quality=True)
        h4 = x.h(p=p4,T=t4,x=x4)

        #Estado 4 real
        h4real=h3+(h4-h3)/nu_p_real

        T4_real,x4_real  = H2O.T_h(h=h4real,p=p4, quality=True)
        s4_real  = H2O.s(T=T4_real,p=p4,x=x4_real)

        ### PARÁMETROS DEL SISTEMA ###

        Wr = h2 - h1 # Trabajo realizado por la turbina
        Wr_real = h2real - h1 #Trabajo real realizado por la turbina
        Q_out = h3-h2 # Calor cedido en el condensador
        Q_out_real = h3-h2real #Calor real cedido por la turbina
        Wp = (h4-h3) #Trabajo de la bomba
        Wp_real = h4real - h3 #Trabajo real de la bomba
        Q_in = h1-h4 #Calor absorbido por el ciclo
        Q_in_real = h1-h4real #Calor real absorbido por el ciclo
        m = Wc/(Q_in_real+Q_out_real) #Flujo másico

        r = abs((Wr-Wp)*100/Q_in) #Rendimiento
        r_real = abs((Wr_real-Wp_real)*100/Q_in_real) #Rendimiento real

        ### REPRESENTACIÓN GRÁFICA ###

        fig, ax = plt.subplots()

        # Del 1 al 2
        T = np.array([t1,t2])
        s = np.array([s1,s2])
        ax.plot(s,T,c='r',label='Ciclo ideal')

        #Del 1 al 2 real
        T = np.array([t1,T2_real])
        s = np.array([s1,s2_real])
        ax.plot(s,T,c='k',ls='--',label='Ciclo real')

        #Del 2 al 3
        ax.plot([s2,s3],[t2,t3],c='r')

        #Del 2 real al 3
        ax.plot([s2_real,s3],[T2_real,t3],c='k',ls='--')

        #del 3 al 4
        T = np.array([t3,t4])
        s = np.array([s3,s4])
        ax.plot(s,T,c='r')

        #del 3 al 4 real
        T = np.array([t3,T4_real])
        s = np.array([s3,s4_real])
        ax.plot(s,T,c='k')

        #del 4 al 1
        T = np.linspace(t4,t1,100)
        p = p4 * np.ones(T.shape)
        s = H2O.s(T=T,p=p)
        s[-1] = s1 # force the last points to be liquid - not vapor
        ax.plot(s,T,c='r')
        ax.set(title='Diagrama T-s Rankine',xlabel='Entropía kJ/kg*K',ylabel='Temperatura (K)')

        #del 4 real al 1
        T = np.linspace(T4_real,t1,100)
        p = p4 * np.ones(T.shape)
        s = H2O.s(T=T,p=p)
        s[-1] = s1 # force the last points to be liquid - not vapor
        ax.plot(s,T,c='k',ls='--')
        ax.set(title='Diagrama T-s Rankine',xlabel='Entropía kJ/kg*K',ylabel='Temperatura (K)')


        # Get the critical and triple point properties
        Tt,pt = H2O.triple()
        Tc,pc = H2O.critical()

        # Explore the temperatures between Tt and Tc in 5K increments
        T = np.linspace(Tt,Tc,100)
        s0 = np.zeros(len(T))
        s_1 = np.zeros(len(T))
        for i in range(len(T)):
            s0[i]=H2O.s(T[i],x=0)
            s_1[i]=H2O.s(T[i],x=1)
        ax.plot(s0,T,'y',s_1,T,'b',ls='--')
        ax.plot(s1,t1,'o',label='Estado 1')
        ax.plot(s2,t2,'o',label='Estado 2')
        ax.plot(s2_real,T2_real,'o',label='Estado 2 Real')
        ax.plot(s3,t3 ,'o',label='Estado 3')
        ax.plot(s4,t4 ,'o',label='Estado 4')
        ax.plot(s4_real,T4_real ,'o',label='Estado 4 Real')
        plt.legend()

        print('RANKINE IDEAL')

        data = pd.DataFrame({'': ['Pressure (MPa)', 'Temperature (ºC)', 'Quality', 'Enthalpy (kJ/kg*K)', 'Entropy (kJ/K)'],
    	'State 1': [p1,t1,x1,h1,s1],
    	'State 2': [p2,t2,x2,h2,s2],
        'State 2r': [p2,T2_real,x2_real,h2real,s2_real],
        'State 3': [p3,t3,x3,h3,s3],
        'State 4': [p4,t4,x4,h4,s4],
        'State 4r': [p4,T4_real,x4_real,h4real,s4_real]})

        print(data.T)

        print('Wp:',Wp)
        print('Wr:',Wr)
        print('Q_in:',Q_in)
        print('Q_out:',Q_out)

        print('Mass flow rate:',Wc/(Q_in+Q_out))
        print('Real mass flow rate',Wc/(Q_in_real+Q_out_real))
        print('Rendimiento {} % | Rendimiento real {} %'.format(r,r_real))
        plt.show()

        return r_real,m

Rankine()

print('--------------------')

#RANKINE SOBRECALENTADO

def Rankine_overheat(x=pm.get('mp.H2O'),p1=p1,p2=p2,x1=x1,x3=0,t1_=450,nu_t_real=0.85,nu_p_real=0.85,w_cycle=200):

    #Estado 1
    h1 = x.h(p=p1,x=x1)
    s1 = x.s(p=p1,x=x1)
    t1 = x.T(p=p1,x=x1)

    #Estado 1 sobrecalentado
    t1_ = t1_

    h1_ = x.h(p=p1,T=t1_)
    s1_ = x.s(p=p1,T=t1_)
    x1_ = x.T_s(s=s1_,p=p1, quality=True)[1]

    #Estado 2 sobrecalentado
    t2, x2 = x.T_s(s=s1,p=p2, quality=True)
    s2_ = s1_

    x2_ = x.T_s(s=s2_,p=p2, quality=True)[1]
    h2_ = x.h(T=t2,x=x2_)

    #Estado 2 sobrecalentado real
    h2real_ = h1_ - nu_t_real*(h1_-h2_)
    T2_real_sobrecalentado,x2_real_sobrecalentado  = H2O.T_h(h=h2real_,p=p2, quality=True)
    s2_real_sobrecalentado  = H2O.s(T=T2_real_sobrecalentado,p=p2,x=x2_real_sobrecalentado)

    #Estado 3
    t3=t2
    x3=x3

    p3 = x.p(T=t3, x=x3)
    s3 = x.s(p=p3, T=t3, x=x3)
    h3 = x.h(p=p3, T=t3, x=x3)

    #Estado 4
    s4=s3
    p4=p1

    t4,x4 = x.T_s(s=s4,p=p4, quality=True)
    h4 = x.h(p=p4,T=t4,x=x4)

    #Estado 4 real
    h4real=h3+(h4-h3)/nu_p_real
    T4_real,x4_real  = H2O.T_h(h=h4real,p=p4, quality=True)
    s4_real  = H2O.s(T=T4_real,p=p4,x=x4_real)

    #PARÁMETROS DEL SISTEMA
    wt = h1_-h2real_ #Trabajo turbinas
    wb = h4real-h3  #Trabajo bomba
    qcald = h1_-h4real #Calor absorbido caldera
    qcond = h2real_-h3 #Calor cedido condensador
    r = (wt-wb)*100/qcald #rendimiento
    m = w_cycle*1000/(qcald-qcond) #Flujo másico

    print('RANKINE SOBRECALENTADO')

    data = pd.DataFrame({'': ['Pressure (MPa)', 'Temperature (ºC)', 'Quality', 'Enthalpy (kJ/kg*K)', 'Entropy (kJ/K)'],
	'State 1': [p1,t1,x1,h1,s1],
	'State 1b': [p1,t1_,x1_,h1_,s1_],
    'State 2': [p2,t2,x2_,h2_,s2_],
    'State 2r': [p2,T2_real_sobrecalentado,x2_real_sobrecalentado,h2real_,s2_real_sobrecalentado],
    'State 3': [p3,t3,x3,h3,s3],
    'State 4': [p4,t4,x4,h4,s4],
    'State 4r': [p4,T4_real,x4_real,h4real,s4_real]})

    print(data.T)

    print('W_tur:',wt)
    print('W_pump:',wb)
    print('Q_cald:',qcald)
    print('Q_cond:',qcond)
    print('Mass flow rate:',m)
    print('Rendimiento {} %'.format(abs(r)))

    #REPRESENTACIÓN GRÁFICA

    #Se podía haber comprimido más como se ha hecho en el ciclo recalentado.

    fig, ax = plt.subplots()

    # Del 1 al 1'
    T = np.array([t1,t1_])
    s = np.array([s1,s1_])
    ax.plot(s,T,c='r',label='Ciclo ideal')

    #Estado 1 al 2'
    T = np.array([t1_,t2])
    s = np.array([s1_,s2_])
    ax.plot(s,T,c='r')

    #Del 1' al 2' real
    T = np.array([t1_,T2_real_sobrecalentado])
    s = np.array([s1_,s2_real_sobrecalentado])
    ax.plot(s,T,c='k',ls='--',label='Ciclo real')

    #del 2' real al 3
    ax.plot([s2_real_sobrecalentado,s3],[T2_real_sobrecalentado,t3],c='k',linewidth=1.5,ls='--')

    # de 2' a 3
    T = np.array([t2,t3])
    s = np.array([s2_,s3])
    ax.plot(s,T,c='r')

    #del 3 al 4
    T = np.array([t3,t4])
    s = np.array([s3,s4])
    ax.plot(s,T,c='r')

    #del 3 al 4 real
    T = np.array([t3,T4_real])
    s = np.array([s3,s4_real])
    ax.plot(s,T,c='k')

    #del 4 al 1
    T = np.linspace(t4,t1,100)
    p = p4 * np.ones(T.shape)
    s = H2O.s(T=T,p=p)
    s[-1] = s1 # force the last points to be liquid - not vapor
    ax.plot(s,T,c='r')
    ax.set(title='Diagrama T-s Rankine',xlabel='Entropía kJ/kg*K',ylabel='Temperatura (K)')

    #del 4 real al 1
    T = np.linspace(T4_real,t1,100)
    p = p4 * np.ones(T.shape)
    s = H2O.s(T=T,p=p)
    s[-1] = s1 # force the last points to be liquid - not vapor
    ax.plot(s,T,c='k',ls='--')
    ax.set(title='Diagrama T-s Rankine sobrecalentado',xlabel='Entropía kJ/kg*K',ylabel='Temperatura (K)')


    # Get the critical and triple point properties
    Tt,pt = H2O.triple()
    Tc,pc = H2O.critical()

    # Explore the temperatures between Tt and Tc in 5K increments
    T = np.linspace(Tt,Tc,100)
    s0 = np.zeros(len(T))
    s_1 = np.zeros(len(T))
    for i in range(len(T)):
        s0[i]=H2O.s(T[i],x=0)
        s_1[i]=H2O.s(T[i],x=1)
    ax.plot(s0,T,'y',s_1,T,'b',ls='--')

    ax.plot(s1,t1,'o',label='Estado 1a')
    ax.plot(s1_,t1_ ,'o',label='Estado 1b')
    ax.plot(s2_,t2,'o',label='Estado 2')
    ax.plot(s2_real_sobrecalentado,T2_real_sobrecalentado,'o',label='Estado 2 Real')
    ax.plot(s3,t3 ,'o',label='Estado 3')
    ax.plot(s4,t4 ,'o',label='Estado 4')
    ax.plot(s4_real,T4_real ,'o',label='Estado 4 Real')

    plt.legend()
    plt.show()

    return r,m

Rankine_overheat()

'''El ciclo sobrecalentado mejora la eficiencia del cilco Rankine ideal.
 Además, el título de vapor aumenata lo cual significa menos cantidad de gotas de agua que puedan dañar la turbina.'''

print('--------------------')

#CICLO RANKINE RECALENTADO

def Rankine_reheat(x=pm.get('mp.H2O'),t1=500, p1 = 10, p2 = 0.8, t3 = 460, p4 = 0.008, x5 = 0, w_cycle = 200, r_t1 = 0.91, r_t2 = 0.93,r_p = 0.88):

    #Estado 1
    t1 = t1
    p1 = p1
    x1 = -1 # x = -1 porque el punto se encuentra por fuera de la campana

    h1 = x.h(p=p1,x=x1,T=t1)
    s1 = x.s(p=p1,x=x1,T=t1)

    #Estado 2
    p2 = p2
    s2 = s1

    t2,x2 = x.T_s(s=s2,p=p2, quality=True)
    h2 = x.h(p=p2,T=t2,x=x2)

    #Estado 2 real
    h2r=h1-(r_t1*(h1-h2))
    T2r,x2r  = H2O.T_h(h=h2r,p=p2, quality=True)
    s2r  = H2O.s(T=T2r,p=p2,x=x2r)

    #Estado 3
    t3 = t3
    p3 = p2

    s3 = x.s(p=p3,T=t3)
    x3 = x.T_s(s=s3,p=p3, quality=True)[1]
    h3 = x.h(p=p3,T=t3,x=x3)

    #Estado 4
    s4 = s3
    p4 = p4

    t4,x4 = x.T_s(s=s4,p=p4,quality=True)
    h4 = x.h(p=p4,T=t4,x=x4)

    #Estado 4 real
    h4r = h3-r_t2*(h3-h4)
    T4r,x4r  = H2O.T_h(h=h4r,p=p4, quality=True)
    s4r  = H2O.s(T=T4r,p=p4,x=x4r)

    #Estado 5
    t5=t4
    x5=x5

    p5 = x.p(T=t5, x=x5)
    s5 = x.s(p=p5, T=t5, x=x5)
    h5 = x.h(p=p5, T=t5, x=x5)

    #Estado 6
    s6=s5
    p6=p1

    t6,x6 = x.T_s(s=s6,p=p6, quality=True)
    h6 = x.h(p=p6,T=t6,x=x6)

    #Estado 6 real
    h6r=h5+(h6-h5)/r_p
    T6r,x6r  = H2O.T_h(h=h6r,p=p6, quality=True)
    s6r  = H2O.s(T=T6r,p=p6,x=x6r)

    ### PARÁMETROS DEL SISTEMA ###

    wt = (h1-h2r)+(h3-h4r) #Trabajo turbinas
    wb = h6r-h5 #Trabajo bomba
    qcald = (h1-h6r)+(h3-h2r) #Calor absorbido caldera
    qcond = h4r-h5 #Calor cedido condensador
    r = (wt-wb)*100/qcald #rendimiento
    m = w_cycle*1000/(qcald-qcond) #Flujo másico

    #Recomendable comentar los print antes de ejecutar el análisis de sensibilidad
    #para que no muestre los parámetros para cada iteración.

    print('RANKINE RECALENTADO')

    data = pd.DataFrame({'': ['Pressure (MPa)', 'Temperature (ºC)', 'Quality', 'Enthalpy (kJ/kg*K)', 'Entropy (kJ/K)'],
	'State 1': [p1,t1,x1,h1,s1],
	'State 2r': [p2,T2r,x2r,h2r,s2r],
    'State 3': [p3,t3,x3,h3,s3],
    'State 4r': [p4,T4r,x4r,h4r,s4r],
    'State 5': [p5,t5,x5,h5,s5],
    'State 6r': [p6,T6r,x6r,h6r,s6r]})

    print(data.T)

    print('Wtur:',wt)
    print('Wpump:',wb)
    print('Q_cald:',qcald)
    print('Q_cond:',qcond)
    print('Mass flow rate:',m)
    print('Rendimiento {} %'.format(abs(r)))

    #Lo dejo comentado para que al ejecutar el análisis de sensibilidad no saque una gráfica para cada estado.

    # #REPRESENTACIÓN GRÁFICA
    #
    #
    #
    # fig, ax = plt.subplots()
    #
    # # Del 1 al 6
    #
    # ax.plot([s1,s2],[t1,t2],[s2,s3],[t2,t3],[s3,s4],[t3,t4],[s4,s5],[t4,t5],[s5,s6],[t5,t6],c='r',label='Ciclo Ideal')
    # ax.plot([s1,s2r],[t1,T2r],[s2r,s3],[T2r,t3],[s3,s4r],[t3,T4r],[s4r,s5],[T4r,t5],[s5,s6r],[t5,T6r],label='Ciclo Real',ls='--',c='k')
    #
    # #del 6 al 1
    #
    # T = np.linspace(T6r,t1,100)
    # p = p1 * np.ones(T.shape)
    # s = x.s(T=T,p=p)
    # ax.plot(s,T,c='r')
    # ax.set(title='Diagrama T-s Rankine',xlabel='Entropía kJ/kg*K',ylabel='Temperatura (K)')
    #
    # # Get the critical and triple point properties
    # Tt,pt = x.triple()
    # Tc,pc = x.critical()
    #
    # # Cruva cambio de fase
    # T = np.linspace(Tt,Tc,100)
    # sx0 = np.zeros(T.shape)
    # sx1 = np.zeros(T.shape)
    # for i in range(len(T)):
    #     sx0[i]=x.s(T[i],x=0)
    #     sx1[i]=x.s(T[i],x=1)
    # ax.plot(sx0,T,'y',sx1,T,'b',ls='--',label='H2O')
    # plt.show()
    #
    # print('wt:',wt)
    # print('wb:',wb)
    # print('qcald:',qcald)
    # print('qcond:',qcond)
    #
    # print('Rendimiento:',r)
    # print('mass flow:',m)

    return r,m

Rankine_reheat()


'''El ciclo recalentado mejora la eficiencia del cilco Rankine sobrecalentado.
 También disminuye el flujo másico para el mismo trabajo realizado del ciclo.
 Además, el título de vapor aumenata lo cual significa menos cantidad de gotas de agua que puedan dañar la turbina.'''



'''
ANÁLISIS DE SENSIBILIDAD DEL CICLO RANKINE RECALENTADO.

Eficiciencia y flujo másico del ciclo en función de las temperaturas T1 y T3.
Es un mapa de calor que muestra como la eficiencicia aumenta con el aumento de ambas temperaturas.
'''

temp_1t = np.linspace(500,700,5) #temperatura de entrada de la primera turbina
temp_2t = np.linspace(460,660,5) #temperatura entrada segunda turbina

def sens_temp(t1,t3):
    r=[]
    m=[]
    t1 = []
    t3 = []
    for i in temp_1t:
        for j in temp_2t:
            r.append(Rankine_reheat(t1 = i,t3 = j)[0][0]) #Aplicable a Rankine() y Rankine_overheat()
            m.append(Rankine_reheat(t1 = i,t3 = j)[1][0]) #Aplicable a Rankine() y Rankine_overheat()
            t1.append(i)
            t3.append(j)

    df = pd.DataFrame({'Temp_entrada_primera_turb':t1,'Temp_entrada_segunda_turbina':t3,'r':r})
    df = df.pivot('Temp_entrada_primera_turb','Temp_entrada_segunda_turbina','r')
    df = np.round(df, 2)


    df1 = pd.DataFrame({'Temp_entrada_primera_turb':t1,'Temp_entrada_segunda_turbina':t3,'m':m})
    df1 = df1.pivot('Temp_entrada_primera_turb','Temp_entrada_segunda_turbina','m')
    df1 = np.round(df1, 2)

    fig,ax = plt.subplots(1,2)
    sns.heatmap(df,annot=True,fmt='.2f',ax=ax[0],cbar_kws={'label': 'Rendimiento (%)'})
    sns.heatmap(df1,annot=True,fmt='.2f',ax=ax[1],cmap="YlGnBu",cbar_kws={'label': 'Flujo másico (kg/s)'})

    return plt.show()

sens_temp(temp_1t,temp_2t)

'''
Eficiciencia y flujo másico del ciclo en función de las presiones P1 y P2.
Mapa de calor que muestra el aumento de la eficiecnia con el aumento de ambas presiones.
Aquí se consigue un disminución de del flujo másico cuanto mayor es la diferencia entre ambas presiones.
'''

pres1 = np.linspace(8,15,5) #presion entrada primera turbina
pres2 = np.linspace(0.8,2,5) #presion salida primera turbina

def sens_pres12(p1,p2):
    r=[]
    m=[]
    p1 = []
    p2 = []
    for i in pres1:
        for j in pres2:
            r.append(Rankine_reheat(p1 = i,p2 = j)[0][0]) #Aplicable a Rankine() y Rankine_overheat()
            m.append(Rankine_reheat(p1 = i,p2 = j)[1][0]) #Aplicable a Rankine() y Rankine_overheat()
            p1.append(i)
            p2.append(j)

    df = pd.DataFrame({'Pres_entrada_primera_turbina':p1,'Expansion_tras_primera_turbina':p2,'r':r})
    df = df.pivot('Pres_entrada_primera_turbina','Expansion_tras_primera_turbina','r')
    df = np.round(df, 2)

    df1 = pd.DataFrame({'Pres_entrada_primera_turbina':p1,'Expansion_tras_primera_turbina':p2,'m':m})
    df1 = df1.pivot('Pres_entrada_primera_turbina','Expansion_tras_primera_turbina','m')
    df1 = np.round(df1, 2)

    fig,ax = plt.subplots(1,2)
    sns.heatmap(df,annot=True,fmt='.2f',ax=ax[0],cbar_kws={'label': 'Rendimiento (%)'})
    sns.heatmap(df1,annot=True,fmt='.2f',ax=ax[1],cmap="YlGnBu",cbar_kws={'label': 'Flujo másico (kg/s)'})

    return plt.show()

sens_pres12(pres1,pres2)

'''
Eficiciencia y flujo másico del ciclo en función de las presiones P1 y P4.
Mapa de calor que muestra que la eficiencia aumenta al disminuir la presión P4,
y también disminuye el flujo másico cuando disminuimos P4.
'''

pres1 = np.linspace(8,15,5) #presion entrada primera turbina de prueba
pres4 = np.linspace(0.008,0.001,5) #presion entrada condensador de prueba

def sens_pres14(p1,p4):
    r=[]
    m=[]
    p1 = []
    p4 = []
    for i in pres1:
        for j in pres4:
            r.append(Rankine_reheat(p1 = i,p4 = j)[0][0]) #Aplicable a Rankine() y Rankine_overheat()
            m.append(Rankine_reheat(p1 = i,p4 = j)[1][0]) #Aplicable a Rankine() y Rankine_overheat()
            p1.append(i)
            p4.append(j)

    df = pd.DataFrame({'Pres_entrada_primera_turbina':p1,'Presion_entrada_condensador':p4,'r':r})
    df = df.pivot('Pres_entrada_primera_turbina','Presion_entrada_condensador','r')
    df = np.round(df, 2)

    df1 = pd.DataFrame({'Pres_entrada_primera_turbina':p1,'Presion_entrada_condensador':p4,'m':m})
    df1 = df1.pivot('Pres_entrada_primera_turbina','Presion_entrada_condensador','m')
    df1 = np.round(df1, 2)

    fig,ax = plt.subplots(1,2)
    sns.heatmap(df,annot=True,fmt='.2f',ax=ax[0],cbar_kws={'label': 'Rendimiento (%)'})
    sns.heatmap(df1,annot=True,fmt='.2f',ax=ax[1],cmap="YlGnBu",cbar_kws={'label': 'Flujo másico (kg/s)'})

    return plt.show()

sens_pres14(pres1,pres4)
