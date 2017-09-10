import numpy as np
import matplotlib.pyplot as plt

class conformism2D() :
    def __init__(self, N, T, A, V, alpha, tstep) :
        self.N              = N
        self.T              = T
        self.phi            = abs(np.random.normal(0, A, size = (self.N, self.N)))
        self.V              = V
        self.alpha          = alpha
        self.tstep          = tstep
        self.teq            = 100
        self.startconfig    = None
        self.config         = None
        self.equilibriumconfig = None
        self.endconfig      = None
        self.energy         = np.zeros(self.tstep + self.teq)
        self.magnetization  = np.zeros(self.tstep + self.teq)
        self.heatcapacity   = np.zeros(self.tstep + self.teq)
        self.susceptibility = np.zeros(self.tstep + self.teq)
        self.acceptance     = np.zeros(self.tstep + self.teq)

    def simulation(self) :
        # Initiate a configuration of spins
        if self.startconfig == None :
            self.config = np.around(0.5*self.alpha + ( ( np.random.randint(2, size = (self.N, self.N)) - 0.5 )*(1+self.alpha) ), decimals = 2)
            self.startconfig = np.copy(self.config)
        else :
            self.config = np.copy(self.startconfig)
        # Run metropolis to reach equilibrium
        for t1 in range(self.teq) :
            self.config, self.acceptance[t1] = self.metropolis(self.config)

        self.equilibriumconfig = np.copy(self.config)
        
        E1 = E2 = M1 = M2 = 0
        
        # Once equilibrium is reached, run metropolis and compute average quantities
        for t2 in range(self.tstep) :
            self.config, self.acceptance[t2+self.teq] = self.metropolis(self.config)
        
        
            if t1%10 == 0 :
                E = self.calcEnergy(self.config)
                M = self.calcMagn(self.config)
                E1 = E1 + E
                E2 = E2 + E*E
                M1 = M1 + M
                M2 = M2 + M*M

                self.energy[self.teq + t2]         = E1 / (self.tstep * self.N*self.N)
                self.magnetization[self.teq + t2]  = M1 / (self.tstep * self.N*self.N)
                self.heatcapacity[self.teq + t2]   = (E2 / self.tstep - E1*E1 / (self.tstep * self.tstep)) / (self.N * self.T * self.T)
                self.susceptibility[self.teq + t2] = (M2 / self.tstep - M1*M1/ (self.tstep * self.tstep)) / (self.N * self.T)
                
        self.endconfig = np.copy(self.config)
            

    def metropolis(self, config) :
        acceptance = 0

        def energy(config, s, i, j) :
            a = float(self.phi[i,j])
            nn = config[(i+1)%self.N, j] + config[(i-1)%self.N, j] + config[i, (j+1)%self.N] + config[i, (j-1)%self.N]
            return (a - 0.5*float(self.V)*nn)*float(s)
        
        for l in range(self.N) :
            for k in range(self.N) :
                i, j = np.random.randint(0, self.N), np.random.randint(0, self.N)
                s = config[i, j]
                e = energy(config, s, i, j)

                if s == -0.5 :
                    stry = 0.6
                else :
                    stry = -0.5
                    
                etry = energy(config, stry, i, j)
                alpha = min(1, np.exp(-(e-etry)/self.T))
                if alpha > np.random.random() :
                    config[i, j] = stry
                    acceptance += 1
                    
        return config, acceptance

    def calcEnergy(self, config) :
        e = 0
        for i in range(self.N) :
            for j in range(self.N) :
                nn = config[(i+1)%self.N, j] + config[(i-1)%self.N, j] + config[i, (j+1)%self.N] + config[i, (j-1)%self.N]
                e += (self.phi[i, j] - 0.5*self.V*nn)*self.config[i, j]
        return e

    def calcMagn(self, config) :
        return np.sum(config)

    
