# %%
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
from qiskit.quantum_info import Pauli

# %% [markdown]
# # Hamiltonia Definitions and Trotterizations

# %% [markdown]
# ### Traverse Field Ising Modelm (TFIM)


class TrotterHamiltonian:
    def __init__(self, T: int, Nk: int, Rq: int, Hk: list, indk: np.ndarray):
        self.T = T          #Trotter size
        self.Nk = Nk        #number of pieces
        self.Rq = Rq        #number of times each single qubit operator acts
        self.Hk = Hk        #list of all the Hammiltonian terms
        self.indk = indk    #first qbit on which each piece acts on

# %%
def TFIM(J,h,n_qubits,T=2,sparse=True):
    """
    This function ...
    """
    X=[]
    ZZ=[]
    M=[]
    for i in range(n_qubits):
        # Hamiltonian
        sx=["I"]*n_qubits
        szz=["I"]*n_qubits
        i1= (i+1)%n_qubits #periodic index
        sx[i]="X"
        szz[i]="Z"
        szz[i1]="Z"
        X.append(Pauli(''.join(sx)).to_matrix(sparse=sparse))
        ZZ.append(Pauli(''.join(szz)).to_matrix(sparse=sparse))
        #Magnetization
        mag=["I"]*n_qubits
        mag[i]="Z"
        M.append(Pauli(''.join(mag)).to_matrix(sparse=sparse))
        
    H = J*np.sum(ZZ,axis=0) + h*np.sum(X,axis=0)
    if sparse:
        H = sp.csc_matrix(H)
    #############################################################
    H_T=[]
    if T>n_qubits:
        print("Your Trotterization size is bigger than you system")
    if T<2:
        print("Your Trotterization size is too small, min value is T=2")
    N_k=Fraction(n_qubits,T).numerator
    R=Fraction(n_qubits,T).denominator
    if R==1:
        N_k=N_k*2
        R=R*2
    index=np.floor(n_qubits/N_k*np.arange(N_k)).astype(int)
    for i in range(N_k):
        ind=index[i]
        h_k=np.zeros((2**n_qubits,2**n_qubits),dtype=complex)
        for j in range(T):
            indT=(ind+j)%n_qubits
            if indT==(index[(i+1)%N_k]-1)%n_qubits:
                h_k=h_k+J*(ZZ[indT])/(R-1)+h*X[indT]/R
                #print(i,indT,"eliminat central")
            elif j == T-1:
                h_k=h_k+h*X[indT]/R
                #print(i,indT,"eliminat final")
            else:
                h_k=h_k+J*(ZZ[indT])/R+h*X[indT]/R
                #print(i,indT,"terme normal")
        if sparse:
            H_T.append(sp.csc_matrix(h_k))
        else:
            H_T.append(h_k)
    
    ###################################################
    #This is a check to see if the trotterization is the same as the original Hamiltonian
    if sparse:
        difH=sp.linalg.norm((H-sum(H_T)),ord=2)
    else:
        difH=np.linalg.norm((H-sum(H_T)),ord=2)
    if difH<1e-14:
        print('Succesfull Troterization')
        print("The Trotterization consists of",N_k,"terms with the starting qubit of each piece at",index)
        print("Each single qubit term appears",R,"times")
    else:
        print('Failed Trotterization, you can still use the generated full Hamiltonian')

    H_trot=TrotterHamiltonian(T,N_k,R,H_T,index)

    return H,H_trot

# %% [markdown]
# ### Glauber Ising Model

# %%
def GI(gam,delta,n_qubits,T=3,Tstr=3,R3=1,sparse=True):
    """
    This function ...
    """
    X=[]
    ZXZ=[]
    ZZ=[]
    ZIZ=[]
    M=[]
    for i in range(n_qubits):
        sixi=[]
        szxz=[]
        szzi=[]
        sizz=[]
        sziz=[]
        mag=[]
        for k in range(n_qubits):
            sixi.append("I")
            szxz.append("I")
            szzi.append("I")
            sizz.append("I")
            sziz.append("I")
            mag.append("I")
        i1= (i+1)%n_qubits
        i2= (i+2)%n_qubits
        sixi[i1]="X"
        
        szxz[i]="Z"
        szxz[i1]="X"
        szxz[i2]="Z"
        
        szzi[i]="Z"
        szzi[i1]="Z"
        
        sizz[i1]="Z"
        sizz[i2]="Z"
        
        sziz[i]="Z"
        sziz[i2]="Z"
        mag[i]="Y"
        xstr=''.join(sixi)
        zxzstr=''.join(szxz)
        zzistr=''.join(szzi)
        zizstr=''.join(sziz)
        izzstr=''.join(sizz)
        magstr=''.join(mag)
        X.append(Pauli(xstr).to_matrix(sparse=sparse))
        ZXZ.append(Pauli(zxzstr).to_matrix(sparse=sparse))
        ZZ.append(Pauli(zzistr).to_matrix(sparse=sparse)+Pauli(izzstr).to_matrix(sparse=sparse))
        #ZZ.append(Pauli(izzstr).to_matrix(sparse=sparse))
        ZIZ.append(Pauli(zizstr).to_matrix(sparse=sparse))
        M.append(Pauli(magstr).to_matrix(sparse=sparse))
    idstr=[]
    for k in range(n_qubits):
        idstr.append("I")
    idstr2=''.join(idstr)
    ID=Pauli(idstr2).to_matrix(sparse=sparse)
    F=Pauli(''.join(["X"]*n_qubits)).to_matrix(sparse=sparse)
    if gam==0:
        A=(1-delta)/2
    else:
        A=(1+delta)*gam**2/(2*(1-np.sqrt(1-gam**2)))-delta
    B= 1-A
    H = -1/2*(A*np.sum(X,axis=0) - B*np.sum(ZXZ,axis=0) + gam/2*(1+delta)*np.sum(ZZ,axis=0) - delta*np.sum(ZIZ,axis=0)-n_qubits*ID)
    if sparse:
        H = sp.csc_matrix(H)

    ###############################################################################################################################
    #Expanation R3 and Tstr
    #R3 represents the number of times 3 body terms appear in your trotterization, most common values are 1, 2 , 4
    #but sometimes (odd lengths and Ts) is site dependent
    #Tstr is the Trotterization strategy, when your chain can be divided into whole pieces (n_qubits/T is a whole number)
    #to include the interaction terms between pieces you have to add orther pieces, Tstr=2 add one of them, Tstr=3 adds two

    H_T=[]
    if T>n_qubits:
        print("Your Trotterization size is bigger than you system")
    if T<3:
        print("Your Trotterization size is too small, min value is T=3")
    N_k=Fraction(n_qubits,T).numerator
    R=Fraction(n_qubits,T).denominator
    if R==1:
        N_k=N_k*Tstr
        R=R*Tstr
    #HD_sin=HD_sin-1/2*(A*X[ind]/T-B*ZXZ[ind]+gamma/2*(1+delta)*(ZZ[ind])-delta*ZIZ[ind]-ID/T)
    index=np.floor(n_qubits/N_k*np.arange(N_k)).astype(int)
    for i in range(N_k):
        ind=index[i]
        h_k=np.zeros((2**n_qubits,2**n_qubits),dtype=complex)
        for j in range(T):
            indT=(ind+j)%n_qubits
            if j == T-1 or j==T-2:
                h_k=h_k-1/2*(A*X[indT] -ID)/R
                #print(i,indT,"eliminat final")
            else:
                h_k=h_k-1/2*(- B*ZXZ[indT] + gam/2*(1+delta)*ZZ[indT] - delta*ZIZ[indT])/R3-1/2*(A*X[indT] -ID)/R
                #print(i,indT,"terme normal")
        if sparse:
            H_T.append(sp.csc_matrix(h_k))
        else:
            H_T.append(h_k)
    
    ###################################################
    #This is a check to see if the trotterization is the same as the original Hamiltonian
    if sparse:
        difH=sp.linalg.norm((H-sum(H_T)),ord=2)
    else:
        difH=np.linalg.norm((H-sum(H_T)),ord=2)
    if difH<1e-14:
        print('Succesfull Troterization')
        print("The Trotterization consists of",N_k,"terms with the starting qubit of each piece at",index)
        print("Each single qubit term appears",R,"times")
    else:
        print('Failed Trotterization, you can still use the generated full Hamiltonian')

    H_trot=TrotterHamiltonian(T,N_k,R,H_T,index)
    
    return H,H_trot