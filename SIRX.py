import numpy as np
from scipy.integrate import ode
from lmfit import minimize, Parameters
import pickle


class SIRXConfirmedModel:

    def __init__(self):
        pass


    # set equation of motion for SIRX dynamics
    def dxdt(self,t,y,eta,rho,kappa,kappa0):
    
        S = y[0]
        I = y[1]
        X = y[2]
        H = y[3]
     
        dy = np.zeros(4)
        dy[0] = -eta*S*I - kappa0*S
        dy[1] = +eta*S*I - rho*I - kappa*I - kappa0*I
        dy[2] = +kappa*I + kappa0*I 
        dy[3] = +kappa0*S +  rho*I
        


        return dy

    def SIRX(self,t, y0, eta, rho, kappa, kappa0, N, I0_factor):
        #provides the numerical solution of SIRX
        X0 = y0 / N
        I0 = X0 * I0_factor
        S0 = 1-X0-I0
        y0 = np.array([S0, I0, X0, 0.0])
        t0 = t[0]

        t = t[1:]

        r = ode(self.dxdt)

        # Runge-Kutta with step size control
        r.set_integrator('dopri5')

        # set initial values
        r.set_initial_value(y0,t0)

        # set transmission rate and recovery rate
        r.set_f_params(eta,rho,kappa,kappa0)

        result = np.zeros((4,len(t)+1))
        result[:,0] = y0

        # loop through all demanded time points
        for it, t_ in enumerate(t):

            # get result of ODE integration
            y = r.integrate(t_)

            # write result to result vector
            result[:,it+1] = y

        return result

    def residual(self, params, x, data):
        #residual the numerical solution minus the values of the real data
        eta = params['eta']
        rho = params['rho']
        kappa = params['kappa']
        kappa0 = params['kappa0']
        I0_factor = params['I0_factor']
        
        N = params['N']

        result = self.SIRX(x, data[0], eta, rho, kappa, kappa0, N, I0_factor)
        X = result[2,:]
        
        residual = X*N - data

        return residual

    def fit(self,t, data,maxfev,params=None,N=None,Nmax=None,method='leastsq'):
        #used to minimize the residual function
        R0 = 6.2
        rho = 1/8
        eta = R0*rho

        if params is None:
            params = Parameters()
            params.add('eta',value=eta,vary=False)
            params.add('rho',value=rho, vary=False)
            params.add('kappa',value=rho, min=0.2)        
            params.add('kappa0',value=rho/2, min=0.)
            params.add('I0_factor',value=10,min=1.0,vary=True)
            varyN = N is None
            if varyN:
                N = 1e7
            if Nmax is None:
                Nmax=115000000
            params.add('N',value=N,min=1000,max=Nmax,vary=varyN)

        if method=='Nelder':
            out = minimize(self.residual, params, args=(t, data, ),
                           method=method,
                        )
        else:
            out = minimize(self.residual, params, args=(t, data, ),
                           maxfev=maxfev,
                           method=method,
                        )
        
        
        return out
    
    
    def fitparam(self, direction):
        #finding the parameters that fit best to our model
        with open(direction, 'rb') as handle:
            mydata = pickle.load(handle)
        
        model = SIRXConfirmedModel()

        paramDict = {i:{} for i in mydata}
        for i in mydata: #loop for countries
            if i != 'Tibet':
                data = mydata[i]['data']
                t = mydata[i]['times']
                while True:
                    if data[0] < 0.2:
                        data = data[1:]
                        t = t[1:]
                    else:
                        break
                Npop = mydata[i]['population']
            
                result = model.fit(t, data, maxfev=1000000,params=None,N=Npop,Nmax=None,method='Nelder')
                
                    
                parameters = dict(result.params.valuesdict())     
                paramDict[i] = parameters
                paramDict[i].update({'X0' : mydata[i]['data'][0]})

                
                paramDict[i].update( {'totDays' : t[-1]} )
    
                
        return paramDict

