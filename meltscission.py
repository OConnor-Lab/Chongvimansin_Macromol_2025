import numpy as np
from netCDF4 import Dataset

class EnergyAnalysis:
    def __init__(self, M, N, ne, Cinf, kB, T, taue, rho, b, Kbend):
        '''
        Initialize the energy analysis class.

        Parameters
        -------------
        M: int
            Number of chains.
        N: int
            Chain length (# of beads).
        ne: int
            Entanglement length (# of beads).
        Cinf: float
            Flory's characteristic ratio.
        kB: float
            Boltzmann constant.
        T: float
            Temperature.
        taue: float
            Entanglement time.
        rho: float
            Monomer number density (M*N/Volume).
        b: float
            Equilibrium bond length.
        Kbend: float
            Bond bending stiffness.
        '''
        self.M = M
        self.N = N
        self.Na = int(M*N)
        self.ne = ne
        self.Cinf = Cinf
        self.kB = kB
        self.T = T
        self.taue = taue
        self.rho = rho
        self.b = b
        self.Kbend = Kbend

    def getVolume(self, DumpPath):
        '''Input the filepath for the dump file'''
        ds = Dataset(DumpPath, 'r')
        a, b, c = ds['cell_lengths'][:].data.T
        alpha, beta, gamma = ds['cell_angles'][:].data.T
        vol = a * b * c * np.sqrt(1 - (np.cos(alpha / 180 * np.pi)) ** 2 -
                                  (np.cos(beta / 180 * np.pi) ** 2) -
                                  (np.cos(gamma / 180 * np.pi) ** 2) +
                                  (2 * np.cos(alpha / 180 * np.pi) *
                                   np.cos(beta / 180 * np.pi) *
                                   np.cos(gamma / 180 * np.pi)))
        vol = np.mean(vol)
        ds.close()
        return vol

    def Wi(self, rate):
        '''Define a function for calculating the Rouse-Wiesenberg number.'''
        erate = np.float(rate)
        wi = (self.taue / self.ne**2) * erate * self.N**2
        return int(wi)

    def ChainEnds(self, fpath, tf):
        '''Index chain ends.'''
        ds = Dataset(fpath + '/dump_draw.nc', 'r')
        be = ds['c_be'][tf, :][:(self.M * self.N)]
        ds.close()
        mol = np.repeat(np.arange(1, self.M+ 1), self.M)
        cee = (np.min(be[np.nonzero(be)]) + np.median(be)) / 2
        inde = np.array(np.where(be < cee))
        ind0 = np.array(np.where(inde == np.where(be == 0)[0]))[1:]
        for ind00 in ind0:
            inde = np.insert(inde, ind00, 0)
            inde[ind00 + 1] = 0
        inde = inde.reshape(int(np.size(inde) / 2), 2)
        return inde

    def BreakCount(self, fpath):
        '''Counting number of breaking events during each strain interval.
           nbreak: int array
                   number of bond broken at each tf
           phin: float array
                 fraction of chains left unbroken at each tf
        '''
        ds = Dataset(fpath + '/dump_draw.nc', 'r')
        t = np.array(ds['time'][:])
        nf = len(t)  # Number of timeframes
        phin = np.zeros(nf)
        nbreak = np.zeros(nf)
        for tf in range(nf):
            inde = self.ChainEnds(fpath, tf)  # index chain ends at tf
            indN = inde[inde[:, 1] - inde[:, 0] == (self.N - 1)]
            if tf == 0:  # At time frame = 0, nothing to compare with
                indeP = inde
            nbreak[tf] += np.size(np.setdiff1d(inde, indeP)) // 2
            indeP = inde
            phin[tf] = len(indN)
        ds.close()
        return nbreak, phin, t
    
    def MwDistribution(self, fpath, fout, fname):
        '''Calculate the chain length of all broken chain at each timeframe.
        nf: int
            number of timeframes
        '''
        ds = Dataset(fpath+'/dump_draw.nc','r')
        nf = len(ds['time'][:].data) # getting number of time frame
        tarray = np.arange(0,nf)
        ds.close()

        Mw_dist = []
        # looping over output time frame
        for tf in tarray:
            CLen = []
            inde = self.ChainEnds(fpath, tf) # chain ends indices
            clen = inde[:,1] - inde[:,0] + 1
            CLen = np.append(CLen, clen[clen!=400])

            Mw_dist += [CLen]

        return nf, Mw_dist
    
    def sscurve(self, fpath, fout, fname):
        # Get extensional stress and strain from thermo output file
        sv = np.loadtxt(fpath + '/sv.out')
        sigEx = -(sv[:,8]-0.5*(sv[:,6]+sv[:,7])), # -sigE = pzz-(pxx+pyy)/2
        eps = sv[:,1] # strain
        np.save(fout + '/eps_' + fname + '.npy') # strain from thermo output (sv.out)
        np.save(fout + '/sigEx_' + fname + '.npy') # stress from thermo output (sv.out)
    
    def getRM(self, fpath):
        '''Get rotation matrix from UEFEX package.'''
        ds = Dataset(fpath + '/dump_draw.nc', 'r')
        t = ds['time'][:]  # time
        nf = np.size(t)
        RM = np.zeros((nf, 3, 3))
        RM[:, 0, 0] = ds['c_rmatrix[1]'][:]
        RM[:, 0, 1] = ds['c_rmatrix[2]'][:]
        RM[:, 0, 2] = ds['c_rmatrix[3]'][:]
        RM[:, 1, 0] = ds['c_rmatrix[4]'][:]
        RM[:, 1, 1] = ds['c_rmatrix[5]'][:]
        RM[:, 1, 2] = ds['c_rmatrix[6]'][:]
        RM[:, 2, 0] = ds['c_rmatrix[7]'][:]
        RM[:, 2, 1] = ds['c_rmatrix[8]'][:]
        RM[:, 2, 2] = ds['c_rmatrix[9]'][:]
        ds.close()
        return RM

    def getframe(self, fpath):
        '''Get flow frames.'''
        ds = Dataset(fpath + '/dump_draw.nc', 'r')
        t = ds['time'][:]  # time
        nf = np.size(t)

        a, b, c = ds['cell_lengths'][:].data.T
        alpha, beta, gamma = ds['cell_angles'][:].data.T
        ds.close()

        lx = a
        xy = b * np.cos(gamma * np.pi / 180.0)
        xz = c * np.cos(beta * np.pi / 180.0)
        ly = (b**2 - xy**2)**0.5
        yz = (b * c * np.cos(alpha * np.pi / 180.0) - xy * xz) / ly
        lz = (c**2 - xz**2 - yz**2)**0.5
        C = np.zeros((nf, 3, 3))

        C[:, 0, 0] = lx
        C[:, 0, 1] = xy
        C[:, 1, 1] = ly
        C[:, 0, 2] = xz
        C[:, 1, 2] = yz
        C[:, 2, 2] = lz

        return C, nf
    
    def unwrap(self, fpath, i):
        '''Unwrap chains (no image flags).'''
        C, _ = self.getframe(fpath)
        ds = Dataset(fpath + '/dump_draw.nc', 'r')
        Xi = ds['coordinates'][i][:(self.M * self.N)]
        ds.close()
        Ci = C[i]
        invCi = np.linalg.inv(Ci)
        Xsi = np.dot(invCi, Xi.reshape(self.Na, 3).T).T.reshape(self.M, self.N, 3)
        dXsi = Xsi[:, 1:] - Xsi[:, :-1]
        n = np.zeros((self.Na, 3)).reshape(self.M, self.N, 3)
        dXsi[dXsi > 0.5] -= 1
        dXsi[dXsi < -0.5] += 1
        Xsui = np.cumsum(dXsi, axis=1)
        Xsui = np.append(np.zeros((self.M, 1, 3)), Xsui, axis=1)
        Xsui += Xsi[:, 0][:, np.newaxis]
        Xui = np.dot(Ci, Xsui.reshape(self.Na, 3).T).T
        return Xui

    def rotate(self, fpath, i):
        '''Rotate flow frame.'''
        C, _ = self.getframe(fpath)
        RM = self.getRM(fpath)
        Ci = C[i]
        Xui = self.unwrap(fpath, i)
        if i == 0:
            # do not rotate at time = 0
            Cfi = Ci
            Xufi = Xui
        else:
            Cfi = np.dot(RM[i], Ci)
            Xufi = np.dot(RM[i].T, Xui.T).T
        return Xufi

    def QE(self, r, K):
        '''Quartic Potential Calculation.'''
        # Define Quartic potential parameters
        Rc = 1.5
        B1 = 0.0
        B2 = -0.7425
        U0 = 92.74467
        epsi = 1.0
        sig = 1.0
        Eq = K * (r - Rc)**2.0 * (r - Rc - B1) * (r - Rc - B2) + U0
        return Eq

    def BondEnergy(self, fpath):
        '''Obtain bonding energy from LAMMPS.'''
        ds = Dataset(fpath + '/dump_draw.nc', 'r')
        pe = ds['c_be'][:][:(self.Na)]  # system's potential energy
        t = ds['time'][:]
        ds.close()
        return pe, t

    def Sdev(self, Energy, vol, t, rate, de, simple=True):
        '''Take derivative of the energy per unit volume per unit strain.'''
        if simple:
            de = de
        else:
            erate = np.float(rate)
            de = np.mean(t[1:] - t[:-1]) * erate
        stress = np.gradient(Energy / vol) / de
        return stress

    def Pade(self, x):
        '''Treloar's approximation of inverse Langevin function.'''
        fx = 3 * x / (1 - 0.6 * x**2 - 0.2 * x**4 - 0.2 * x**6)
        return fx

    def Theta(self, dXn, aid):
        '''Calculate the bond angles.'''
        norms_aid = np.sqrt(np.einsum('ij,ij->i', dXn[aid], dXn[aid]))
        norms_aid_plus1 = np.sqrt(np.einsum('ij,ij->i', dXn[aid + 1], dXn[aid + 1]))
        rn = norms_aid * norms_aid_plus1
        rdot = np.einsum('ij,ij->i', dXn[aid], dXn[aid + 1])
        costht = np.divide(rdot, rn)
        tht = np.arccos(costht)
        return tht

    def EnergyParameters(self, fpath, fout, fname):
        '''Calculate the parameters needed for energy analysis. Save the parameters to fout as a .np file.'''
        _, nf = self.getframe(fpath) # number of output frames
        ds = Dataset(fpath + '/dump_draw.nc', 'r')
        Ub = ds['c_be'][:].data  # bond PE energy
        ds.close()

        Tht = np.zeros(nf)
        Blen = np.zeros(nf)
        PEquart = np.zeros(nf)
        PEangle = np.zeros(nf)
        sigS = np.zeros(nf)
        ExtRatio = np.zeros(nf)

        for tf in range(nf):
            Xufi = self.rotate(fpath,i=tf)  # rotated coordinations
            inde = self.ChainEnds(fpath, tf)  # chain ends indices
            indN = inde[inde[:, 1] - inde[:, 0] == (self.N - 1)] 
            # only the unbroken chains
            XN = Xufi.reshape(self.M, self.N, 3)[(indN[:, 0] / self.N).astype(int), :, :]
            # at Ne scale
            XN = XN[:, :(self.N // self.ne) * self.ne, :].reshape(indN.shape[0], self.N // self.ne, self.ne, 3)
            dXN = XN[:, :, 1:, :] - XN[:, :, :-1, :] # bond vectors
            dxn = dXN / np.linalg.norm(dXN, axis=3)[..., np.newaxis] # bond unit vectors
            # Calculate orientation tensor
            Otensor = np.einsum('ijkl,ijkm->ijklm', dxn[:, :, :-1], dxn[:, :, 1:])
            Otensor_Ne = np.mean(Otensor, axis=2)

            RN = np.sum(dXN, axis=2)
            rN = np.linalg.norm(RN, axis=2)
            cosN = np.divide(RN[:, :, 2], rN) # Cosine function of the angle between 2 bond vectors
            Lcon = np.sum(np.linalg.norm(dXN, axis=3), axis=2) # contour length
            hN = rN / Lcon # extension ratio
            invL = self.Pade(hN) # inverse Langevin function of hN
            fN = self.kB * self.T / (self.Cinf * self.b) * invL # Langevin force

            # Calculate the Entropic stress
            sig_term = self.rho * hN * fN
            sigNe = sig_term[..., np.newaxis, np.newaxis] * Otensor_Ne
            sigSi = sigNe[:, :, 2, 2] - 0.5 * (sigNe[:, :, 0, 0] + sigNe[:, :, 1, 1])

            # Calculate ...
            bond_vect = Xufi[1:, :] - Xufi[:-1, :]
            ends = np.delete(np.arange(self.Na), inde[:, 1]) 
            iBlen = np.linalg.norm(bond_vect[ends, :], axis=1)

            # Calculate the bond angle potential energy
            angle_id = np.delete(np.arange(self.Na), np.append(inde[:, 1], inde[:, 1] - 1))
            iTht = self.Theta(bond_vect, angle_id)
            iPEangle = self.Kbend * (1 - np.cos(iTht))

            Ub[tf, inde] = 0
            sigS[tf] = np.mean(sigSi)
            PEangle[tf] = np.sum(iPEangle)
            PEquart[tf] = np.mean(self.QE(iBlen, 1500))
            Tht[tf] = np.mean(iTht)
            Blen[tf] = np.mean(iBlen)
            ExtRatio[tf] = np.mean(hN)

        Ub_arr = np.sum(Ub, axis=1) # total Ub

        np.save(fout + '/sigS_ne_' + fname, sigS)
        np.save(fout + '/hN_ne_' + fname, ExtRatio)
        np.save(fout + '/Tht_' + fname, Tht)
        np.save(fout + '/Bl_' + fname, Blen)
        np.save(fout + '/Uq_' + fname, PEquart)
        np.save(fout + '/Ua_' + fname, PEangle)
        np.save(fout + '/Ub_' + fname, Ub_arr)